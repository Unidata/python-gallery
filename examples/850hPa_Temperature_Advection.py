"""
=============================
850 hPa Temperature Advection
=============================

Plot an 850 hPa map with calculating advection using MetPy.

Beyond just plotting 850-hPa level data, this uses calculations from `metpy.calc` to find
the temperature advection. Currently, this needs an extra helper function to calculate
the distance between lat/lon grid points.
"""

###############################
# Imports
from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from metpy.calc import advection
from metpy.units import units
from netCDF4 import num2date
import numpy as np
import scipy.ndimage as ndimage
from siphon.ncss import NCSS


#################################
# Helper functions

# Helper function for finding proper time variable
def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)


# Helper function to calculate distance between lat/lon points
# to be used in differencing calculations
def calc_dx_dy(longitude, latitude, shape='sphere', radius=6370997.):
    """ This definition calculates the distance between grid points that are in
        a latitude/longitude format.

        Using pyproj GEOD; different Earth Shapes
        https://jswhit.github.io/pyproj/pyproj.Geod-class.html

        Common shapes: 'sphere', 'WGS84', 'GRS80'

        Accepts, 1D or 2D arrays for latitude and longitude

        Assumes [Y, X] for 2D arrays

        Returns: dx, dy; 2D arrays of distances between grid points
                 in the x and y direction with units of meters and sign of differencing
    """
    import numpy as np
    from metpy.units import units
    from pyproj import Geod

    if radius != 6370997.:
        g = Geod(a=radius, b=radius)
    else:
        g = Geod(ellps=shape)

    if latitude.ndim == 1:
        longitude, latitude = np.meshgrid(longitude, latitude)

    dy = np.zeros(latitude.shape)
    dx = np.zeros(longitude.shape)

    for i in range(longitude.shape[1]):
        for j in range(latitude.shape[0]-1):
            _, _, dy[j, i] = g.inv(longitude[j, i], latitude[j, i],
                                   longitude[j+1, i], latitude[j+1, i])
    dy[j+1, :] = dy[j, :]

    for i in range(longitude.shape[1]-1):
        for j in range(latitude.shape[0]):
            _, _, dx[j, i] = g.inv(longitude[j, i], latitude[j, i],
                                   longitude[j, i+1], latitude[j, i+1])
    dx[:, i+1] = dx[:, i]

    xdiff_sign = np.sign(longitude[0, 1] - longitude[0, 0])
    ydiff_sign = np.sign(latitude[1, 0] - latitude[0, 0])
    return xdiff_sign*dx*units.meter, ydiff_sign*dy*units.meter


###############################################
# Create NCSS object to access the NetcdfSubset
# ---------------------------------------------
# Data from NOMADS GFS 0.5 deg Analysis Archive
# https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs
dt = datetime(2017, 4, 5, 12)
ncss = NCSS('https://nomads.ncdc.noaa.gov/thredds/ncss/grid/gfs-004-anl/'
            '{0:%Y%m}/{0:%Y%m%d}/gfsanl_4_{0:%Y%m%d}_{0:%H}00_000.grb2'.format(dt))

# Create lat/lon box for location you want to get data for
query = ncss.query().time(dt)
query.lonlat_box(north=65, south=15, east=310, west=220)
query.accept('netcdf4')

# Request data for vorticity
query.variables('Geopotential_height', 'Temperature',
                'U-component_of_wind', 'V-component_of_wind')
data = ncss.get_data(query)

# Pull out variables you want to use
hght_var = data.variables['Geopotential_height']
temp_var = data.variables['Temperature']
u_wind_var = data.variables['U-component_of_wind']
v_wind_var = data.variables['V-component_of_wind']
time_var = data.variables[find_time_var(temp_var)]
lat_var = data.variables['lat']
lon_var = data.variables['lon']

# Get actual data values and remove any size 1 dimensions
lat = lat_var[:].squeeze()
lon = lon_var[:].squeeze()
hght = hght_var[:].squeeze()
temp = temp_var[:].squeeze() * units.K
u_wind = u_wind_var[:].squeeze() * units('m/s')
v_wind = v_wind_var[:].squeeze() * units('m/s')

# Convert number of hours since the reference time into an actual date
time = num2date(time_var[:].squeeze(), time_var.units)

lev_850 = np.where(data.variables['pressure'][:] == 850*100)[0][0]
hght_850 = hght[lev_850]
temp_850 = temp[lev_850]
u_wind_850 = u_wind[lev_850]
v_wind_850 = v_wind[lev_850]

# Combine 1D latitude and longitudes into a 2D grid of locations
lon_2d, lat_2d = np.meshgrid(lon, lat)
# Gridshift for barbs
lon_2d[lon_2d > 180] = lon_2d[lon_2d > 180] - 360


###############################################
# Begin data calculations
# -----------------------

# Use helper function defined above to calculate distance
# between lat/lon grid points
dx, dy = calc_dx_dy(lon_var, lat_var)

# Calculate temperature advection using metpy function
adv = advection(temp_850.T * units.kelvin, [u_wind_850.T, v_wind_850.T],
                (dx.T, dy.T)).T * units('K/sec')

# Smooth heights and advection a little
# Be sure to only put in a 2D lat/lon or Y/X array for smoothing
Z_850 = ndimage.gaussian_filter(hght_850, sigma=3, order=0) * units.meter
adv = ndimage.gaussian_filter(adv, sigma=3, order=0) * units('K/sec')

###############################################
# Begin map creation
# ------------------

# Set Projection of Data
datacrs = ccrs.PlateCarree()

# Set Projection of Plot
plotcrs = ccrs.LambertConformal(central_latitude=[30, 60], central_longitude=-100)

# Create new figure
fig = plt.figure(figsize=(11, 8.5))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                       hspace=0.01, wspace=0.01)

# Add the map and set the extent
ax = plt.subplot(gs[0], projection=plotcrs)
plt.title('850mb Temperature Advection for {0:%d %B %Y %H:%MZ}'.format(time), fontsize=16)
ax.set_extent([235., 290., 20., 55.])

# Add state boundaries to plot
states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                name='admin_1_states_provinces_lakes',
                                                scale='50m', facecolor='none')
ax.add_feature(states_provinces, edgecolor='black', linewidth=1)

# Add country borders to plot
country_borders = cfeature.NaturalEarthFeature(category='cultural',
                                               name='admin_0_countries',
                                               scale='50m', facecolor='none')
ax.add_feature(country_borders, edgecolor='black', linewidth=1)

# Plot Height Contours
clev850 = np.arange(900, 3000, 30)
cs = ax.contour(lon_2d, lat_2d, Z_850, clev850, colors='black', linewidths=1.5,
                linestyles='solid', transform=datacrs)
plt.clabel(cs, fontsize=10, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# Plot Temperature Contours
clevtemp850 = np.arange(-20, 20, 2)
cs2 = ax.contour(lon_2d, lat_2d, temp_850.to(units('degC')), clevtemp850,
                 colors='grey', linewidths=1.25, linestyles='dashed',
                 transform=datacrs)
plt.clabel(cs2, fontsize=10, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

# Plot Colorfill of Temperature Advection
cint = np.arange(-8, 9)
cf = ax.contourf(lon_2d, lat_2d, 3*adv.to(units('delta_degC/hour')), cint[cint != 0],
                 extend='both', cmap='bwr', transform=datacrs)
cax = plt.subplot(gs[1])
cb = plt.colorbar(cf, cax=cax, orientation='horizontal', extendrect=True, ticks=cint)
cb.set_label(r'$^{o}C/3h$', size='large')

# Plot Wind Barbs
ax.barbs(lon_2d, lat_2d, u_wind_850.magnitude, v_wind_850.magnitude,
         length=6, regrid_shape=20, pivot='middle', transform=datacrs)

gs.tight_layout(fig)
plt.show()
