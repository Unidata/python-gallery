"""
======================
Miller Composite Chart
======================

Create a Miller Composite chart based on Miller 1972 in Python with MetPy and
Matplotlib.
"""
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.lines as lines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import metpy.calc as mcalc
from metpy.units import units
from netCDF4 import num2date
import numpy as np
import numpy.ma as ma
from scipy.ndimage import gaussian_filter
from siphon.ncss import NCSS

###########################
# **Get the data**
#
# This example will use data from the North American Mesoscale Model Analysis
# (https://nomads.ncdc.gov/) for 12 UTC 27 April 2011.
ncss = NCSS('https://nomads.ncdc.noaa.gov/thredds/ncss/grid/namanl/201104/20110427/'
            'namanl_218_20110427_1800_000.grb')

# Query for required variables
gfsdata = ncss.query().all_times()
gfsdata.variables('Geopotential_height',
                  'u_wind',
                  'v_wind',
                  'Temperature',
                  'Relative_humidity',
                  'Best_4-layer_lifted_index',
                  'Absolute_vorticity',
                  'Pressure_reduced_to_MSL',
                  'Dew_point_temperature'
                  ).add_lonlat()

# Set the lat/lon box for the data to pull in.
gfsdata.lonlat_box(-135, -60, 15, 65)

# Actually getting the data
data = ncss.get_data(gfsdata)

# Assign variable names to collected data
dtime = data.variables['Geopotential_height'].dimensions[0]
dlev = data.variables['Geopotential_height'].dimensions[1]
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
lev = data.variables[dlev][:] * units.hPa
times = data.variables[dtime]
vtimes = num2date(times[:], times.units)
temps = data.variables['Temperature']
tmp = temps[0, :] * units.kelvin
uwnd = data.variables['u_wind'][0, :] * units.meter / units.second
vwnd = data.variables['v_wind'][0, :] * units.meter / units.second
hgt = data.variables['Geopotential_height'][0, :] * units.meter
relh = data.variables['Relative_humidity'][0, :]
lifted_index = (data.variables['Best_4-layer_lifted_index'][0, 0, :] *
                units(data.variables['Best_4-layer_lifted_index'].units))
Td_sfc = (data.variables['Dew_point_temperature'][0, 0, :] *
          units(data.variables['Dew_point_temperature'].units))
avor = data.variables['Absolute_vorticity'][0, :] * units('1/s')
pmsl = (data.variables['Pressure_reduced_to_MSL'][0, :] *
        units(data.variables['Pressure_reduced_to_MSL'].units))

########################
# Query for 00 UTC to calculate pressure falls and height change

ncss2 = NCSS('https://nomads.ncdc.noaa.gov/thredds/ncss/grid/namanl/201104/20110427/'
             'namanl_218_20110427_0600_000.grb')

# Query for required variables
gfsdata2 = ncss2.query().all_times()
gfsdata2.variables('Geopotential_height',
                   'Pressure_reduced_to_MSL').add_lonlat()

# Set the lat/lon box for the data you want to pull in.
gfsdata2.lonlat_box(-135, -60, 15, 65)

# Actually getting the data
data2 = ncss2.get_data(gfsdata)

hgt_00z = data2.variables['Geopotential_height'][0, :] * units.meter
pmsl_00z = (data2.variables['Pressure_reduced_to_MSL'][0, :] *
            units(data2.variables['Pressure_reduced_to_MSL'].units))

###########################
# **Subset the Data**
#
# With the data pulled in, we will now subset to the specific levels desired

# 300 hPa, index 28
u_300 = uwnd[28, :].to('kt')
v_300 = vwnd[28, :].to('kt')

# 500 hPa, index 20
avor_500 = avor[3, :]
u_500 = uwnd[20, :].to('kt')
v_500 = vwnd[20, :].to('kt')
hgt_500 = hgt[20, :]
hgt_500_00z = hgt_00z[20, :]

# 700 hPa, index 12
tmp_700 = tmp[12, :].to('degC')
rh_700 = relh[12, :]
u_700 = uwnd[12, :].to('kt')
v_700 = vwnd[12, :].to('kt')

# 850 hPa, index 6
tmp_850 = tmp[6, :].to('degC')
u_850 = uwnd[6, :].to('kt')
v_850 = vwnd[6, :].to('kt')
rh_850 = relh[6, :]

##############################
# Calculation of advections will require the help of the following function to find deltas


def calc_dx_dy(longitude, latitude, shape='sphere', radius=6370997.):
    """ This definition calculates the distance between grid points that are in
        a latitude/longitude format.

        Using pyproj GEOD; different Earth Shapes
        https://jswhit.github.io/pyproj/pyproj.Geod-class.html

        Common shapes: 'sphere', 'WGS84', 'GRS80'

        Accepts, 1D or 2D arrays for latitude and longitude

        Assumes [Y, X] for 2D arrays

        Returns: dx, dy; 2D arrays of distances between grid points
                 in the x and y direction with units of meters
    """
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

    xdiff_sign = np.sign(longitude[0, 1]-longitude[0, 0])
    ydiff_sign = np.sign(latitude[1, 0]-latitude[0, 0])
    return xdiff_sign*dx*units.meter, ydiff_sign*dy*units.meter

########################################
# **Prepare Variables for Plotting**
#
# With the data queried and subset, we will make any needed calculations in preparation for
# plotting.
#
# The following fields should be plotted:
#   500-hPa cyclonic vorticity advection
#
#   Surface-based Lifted Index
#
#   The axis of the 300-hPa, 500-hPa, and 850-hPa jets
#
#   Surface dewpoint
#
#   700-hPa dewpoint depression
#
#   12-hr surface pressure falls and 500-hPa height changes


# 500 hPa CVA
dx, dy = calc_dx_dy(lon, lat)
vort_adv_500 = mcalc.advection(avor_500, [v_500.to('m/s'), u_500.to('m/s')], (dy, dx),
                               dim_order='yx') * 1e9
vort_adv_500_smooth = gaussian_filter(vort_adv_500, 4)

####################################
# For the jet axes, we will calculate the windspeed at each level, and plot the highest values
wspd_300 = gaussian_filter(mcalc.get_wind_speed(u_300, v_300), 5)
wspd_500 = gaussian_filter(mcalc.get_wind_speed(u_500, v_500), 5)
wspd_850 = gaussian_filter(mcalc.get_wind_speed(u_850, v_850), 5)

#################################
# 850-hPa dewpoint will be calculated from RH and temperature
Td_850 = mcalc.dewpoint_rh(tmp_850, rh_850 / 100.)

################################
# 700-hPa dewpoint depression will be calculated from temperature and RH
Td_dep_700 = tmp_700 - mcalc.dewpoint_rh(tmp_700, rh_700 / 100.)

######################################
# 12-hr surface pressure falls and 500-hPa height changes
pmsl_change = pmsl - pmsl_00z
hgt_500_change = hgt_500 - hgt_500_00z

######################################
# To plot the jet axes, we will mask the wind fields below the upper 1/3 of windspeed.

mask_500 = ma.masked_less_equal(wspd_500, 0.66 * np.max(wspd_500)).mask
u_500[mask_500] = np.nan
v_500[mask_500] = np.nan

# 300 hPa
mask_300 = ma.masked_less_equal(wspd_300, 0.66 * np.max(wspd_300)).mask
u_300[mask_300] = np.nan
v_300[mask_300] = np.nan

# 850 hPa
mask_850 = ma.masked_less_equal(wspd_850, 0.66 * np.max(wspd_850)).mask
u_850[mask_850] = np.nan
v_850[mask_850] = np.nan

################################
# **Create the Plot**
#
# With the data now ready, we will create the plot


# Set up our projection
crs = ccrs.LambertConformal(central_longitude=-100.0, central_latitude=45.0)


# Coordinates to limit map area
bounds = [(-122., -75., 25., 50.)]
# Choose a level to plot, in this case 296 K
level = 0

# Get data to plot state and province boundaries
states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                name='admin_1_states_provinces_lakes',
                                                scale='50m',
                                                facecolor='none')

#########################
# Plot the composite
fig = plt.figure(1, figsize=(17., 12.))
ax = plt.subplot(111, projection=crs)
ax.set_extent(*bounds, crs=ccrs.PlateCarree())
ax.coastlines('50m', edgecolor='black', linewidth=0.75)
ax.add_feature(states_provinces, edgecolor='black', linewidth=0.25)

# Plot Lifted Index
cs1 = ax.contour(lon, lat, lifted_index, range(-8, -2, 2), transform=ccrs.PlateCarree(),
                 colors='red', linewidths=0.75, linestyles='solid', zorder=7,
                 label='Best Lifted Index')
plt.clabel(cs1, fontsize=10, inline=1, inline_spacing=7,
           fmt='%i', rightside_up=True, use_clabeltext=True)

# Plot Surface pressure falls
cs2 = ax.contour(lon, lat, pmsl_change.to('hPa'), range(-10, -1, 4),
                 transform=ccrs.PlateCarree(),
                 colors='k', linewidths=0.75, linestyles='dashed', zorder=6)
plt.clabel(cs2, fontsize=10, inline=1, inline_spacing=7,
           fmt='%i', rightside_up=True, use_clabeltext=True)

# Plot 500-hPa height falls
cs3 = ax.contour(lon, lat, hgt_500_change, range(-60, -29, 15),
                 colors='k', linewidths=0.75, linestyles='solid', zorder=5)
plt.clabel(cs3, fontsize=10, inline=1, inline_spacing=7,
           fmt='%i', rightside_up=True, use_clabeltext=True)

# Plot surface pressure
ax.contourf(lon, lat, pmsl.to('hPa'), range(990, 1011, 20), alpha=0.5,
            transform=ccrs.PlateCarree(),
            colors='yellow', zorder=1)

# Plot surface dewpoint
ax.contourf(lon, lat, Td_sfc.to('degF'), range(65, 76, 10), alpha=0.4,
            transform=ccrs.PlateCarree(),
            colors=['green'], zorder=2)

# Plot 700-hPa dewpoint depression
ax.contourf(lon, lat, Td_dep_700, range(15, 46, 30), alpha=0.5, transform=ccrs.PlateCarree(),
            colors='tan', zorder=3)

# Plot Vorticity Advection
ax.contourf(lon, lat, vort_adv_500_smooth, range(5, 106, 100), alpha=0.5,
            transform=ccrs.PlateCarree(),
            colors='BlueViolet', zorder=4)

# 300-hPa wind barbs
jet300 = ax.barbs(lon, lat, u_300.m, v_300.m, length=6, regrid_shape=20,
                  transform=ccrs.PlateCarree(),
                  color='green', zorder=10, label='300-hPa Jet Core Winds (kt)')


# 500-hPa wind barbs
jet500 = ax.barbs(lon, lat, u_500.m, v_500.m, length=6, regrid_shape=20,
                  transform=ccrs.PlateCarree(),
                  color='blue', zorder=9, label='500-hPa Jet Core Winds (kt)')

# 850-hPa wind barbs
jet850 = ax.barbs(lon, lat, u_850.m, v_850.m, length=6, regrid_shape=20,
                  transform=ccrs.PlateCarree(),
                  color='k', zorder=8, label='850-hPa Jet Core Winds (kt)')

# Legend
purple = mpatches.Patch(color='BlueViolet', label='Cyclonic Absolute Vorticity Advection')
yellow = mpatches.Patch(color='yellow', label='Surface MSLP < 1010 hPa')
green = mpatches.Patch(color='green', label='Surface Td > 65 F')
tan = mpatches.Patch(color='tan', label='700 hPa Dewpoint Depression > 15 C')
red_line = lines.Line2D([], [], color='red', label='Best Lifted Index (C)')
dashed_black_line = lines.Line2D([], [], linestyle='dashed', color='k',
                                 label='12-hr Surface Pressure Falls (hPa)')
black_line = lines.Line2D([], [], linestyle='solid', color='k',
                          label='12-hr 500-hPa Height Falls (m)')
plt.legend(handles=[jet300, jet500, jet850, dashed_black_line, black_line, red_line,
                    purple, tan, green, yellow], loc=3,
           title='Composite Analysis Valid: {:s}'.format(str(vtimes[0])))

plt.show()
