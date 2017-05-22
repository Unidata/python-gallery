"""
=================================
Geostrophic and Ageostrophic Wind
=================================

Plot a 1000-hPa map calculating the geostrophic from MetPy and finding the
ageostrophic wind from the total wind and the geostrophic wind.

This uses the geostrophic wind calculation from `metpy.calc` to find
the geostrophic wind, then performs the simple subtraction to find the ageostrophic
wind. Currently, this needs an extra helper function to calculate
the distance between lat/lon grid points.

Additionally, we utilize the `ndimage.zoom` method for smoothing the 1000-hPa
height contours without smoothing the data.
"""

########################################
# Imports
from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from metpy.calc import coriolis_parameter, geostrophic_wind
from metpy.units import units
from netCDF4 import num2date
import numpy as np
import scipy.ndimage as ndimage
from siphon.ncss import NCSS

########################################
# Helper function to calculate distance between lat/lon
# -----------------------------------------------------


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

######################################
# Set up access to the data


# Create NCSS object to access the NetcdfSubset
ncss = NCSS('https://nomads.ncdc.noaa.gov/thredds/ncss/grid/gfs-004-anl/201608/20160822/'
            'gfsanl_4_20160822_1800_003.grb2')

# Create lat/lon box for location you want to get data for
query = ncss.query()
query.lonlat_box(north=50, south=30, east=-80, west=-115)
query.time(datetime(2016, 8, 22, 21))

# Request data for geopotential height
query.variables('Geopotential_height', 'U-component_of_wind', 'V-component_of_wind')
query.vertical_level(100000)
data = ncss.get_data(query)

# Pull out variables you want to use
height_var = data.variables['Geopotential_height']
u_wind_var = data.variables['U-component_of_wind']
v_wind_var = data.variables['V-component_of_wind']

# Find the name of the time coordinate
for name in height_var.coordinates.split():
    if name.startswith('time'):
        time_var_name = name
        break
time_var = data.variables[time_var_name]

lat_var = data.variables['lat']
lon_var = data.variables['lon']

# Get actual data values and remove any size 1 dimensions
lat = lat_var[:].squeeze()
lon = lon_var[:].squeeze()
height = height_var[0, 0, :, :].squeeze()
u_wind = u_wind_var[0, 0, :, :].squeeze() * units('m/s')
v_wind = v_wind_var[0, 0, :, :].squeeze() * units('m/s')

# Convert number of hours since the reference time into an actual date
time = num2date(time_var[:].squeeze(), time_var.units)

# Combine 1D latitude and longitudes into a 2D grid of locations
lon_2d, lat_2d = np.meshgrid(lon, lat)

# Smooth height data
height = ndimage.gaussian_filter(height, sigma=1.5, order=0)

# Set up some constants based on our projection, including the Coriolis parameter and
# grid spacing, converting lon/lat spacing to Cartesian
f = coriolis_parameter(np.deg2rad(lat_2d)).to('1/s')
dx, dy = calc_dx_dy(lon_2d, lat_2d)

# In MetPy 0.5, geostrophic_wind() assumes the order of the dimensions is (X, Y),
# so we need to transpose from the input data, which are ordered lat (y), lon (x).
# Once we get the components,transpose again so they match our original data.
geo_wind_u, geo_wind_v = geostrophic_wind(height.T * units.m, f.T, dx.T, dy.T)
geo_wind_u = geo_wind_u.T
geo_wind_v = geo_wind_v.T

# Calculate ageostrophic wind components
ageo_wind_u = u_wind - geo_wind_u
ageo_wind_v = v_wind - geo_wind_v

# Create new figure
fig = plt.figure(figsize=(15, 10), facecolor='black')

# Add the map and set the extent
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-105., -93., 35., 43.])
ax.background_patch.set_fill(False)

# Add state boundaries to plot
states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                name='admin_1_states_provinces_lines',
                                                scale='50m', facecolor='none')
ax.add_feature(states_provinces, edgecolor='white', linewidth=2)

# Add country borders to plot
country_borders = cfeature.NaturalEarthFeature(category='cultural',
                                               name='admin_0_countries',
                                               scale='50m', facecolor='none')
ax.add_feature(country_borders, edgecolor='white', linewidth=2)

# Contour the heights every 10 m
contours = np.arange(10, 200, 10)
# Because we have a very local graphics area, the contours have joints
# to smooth those out we can use `ndimage.zoom`
zoom_500 = ndimage.zoom(height, 5)
zlon = ndimage.zoom(lon_2d, 5)
zlat = ndimage.zoom(lat_2d, 5)
c = ax.contour(zlon, zlat, zoom_500, levels=contours,
               colors='red', linewidth=4)
ax.clabel(c, fontsize=12, inline=1, inline_spacing=3, fmt='%i')

# Set up parameters for quiver plot. The slices below are used to subset the data (here
# taking every 4th point in x and y). The quiver_kwargs are parameters to control the
# appearance of the quiver so that they stay consistent between the calls.
quiver_slices = (slice(None, None, 2), slice(None, None, 2))
quiver_kwargs = {'headlength': 4, 'headwidth': 3, 'angles': 'uv', 'scale_units': 'xy',
                 'scale': 20}

# Plot the wind vectors
wind = ax.quiver(lon_2d[quiver_slices], lat_2d[quiver_slices],
                 u_wind[quiver_slices], v_wind[quiver_slices],
                 color='blue', **quiver_kwargs)
geo = ax.quiver(lon_2d[quiver_slices], lat_2d[quiver_slices],
                geo_wind_u[quiver_slices], geo_wind_v[quiver_slices],
                color='darkorchid', **quiver_kwargs)
ageo = ax.quiver(lon_2d[quiver_slices], lat_2d[quiver_slices],
                 ageo_wind_u[quiver_slices], ageo_wind_v[quiver_slices],
                 color='lime', **quiver_kwargs)

# Add a title to the plot
plt.title('1000mb Geopotential Heights(m), Wind(blue), Geostrophic Wind(purple), and \n'
          'Ageostrophic Wind(green) for {0:%d %B %Y %H:%MZ}'.format(time),
          color='white', size=14)
plt.tight_layout()
plt.show()
