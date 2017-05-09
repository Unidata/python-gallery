"""
==================
Smoothing Contours
==================

Demonstrate how to smooth contour values from a higher resolution
model field.

By: Kevin Goebbert

Date: 13 April 2017
"""

##################################
# Do the needed imports
from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from metpy.calc import get_wind_speed
from metpy.units import units
from netCDF4 import num2date
import numpy as np
import scipy.ndimage as ndimage
from siphon.ncss import NCSS


##################################
# Set up netCDF Subset Service link
dt = datetime(2016, 4, 16, 18)
ncss = NCSS('http://nomads.ncdc.noaa.gov/thredds/ncss/grid/namanl/'
            '{0:%Y%m}/{0:%Y%m%d}/namanl_218_{0:%Y%m%d}_{0:%H}00_000.grb'.format(dt))

# Data Query
hgt = ncss.query().time(dt)
hgt.variables('Geopotential_height', 'u_wind', 'v_wind').add_lonlat()

# Actually getting the data
data = ncss.get_data(hgt)


##################################
# Pull apart the data

# Get dimension names to pull appropriate variables
dtime = data.variables['Geopotential_height'].dimensions[0]
dlev = data.variables['Geopotential_height'].dimensions[1]
dlat = data.variables['Geopotential_height'].dimensions[2]
dlon = data.variables['Geopotential_height'].dimensions[3]

# Get lat and lon data, as well as time data and metadata
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
lons[lons > 180] = lons[lons > 180] - 360

# Need 2D lat/lons for plotting, do so if necessary
if lats.ndim < 2:
    lons, lats = np.meshgrid(lons, lats)

# Determine the level of 500 hPa
levs = data.variables[dlev][:]
lev_500 = np.where(levs == 500)[0][0]

# Create more useable times for output
times = data.variables[dtime]
vtimes = num2date(times[:], times.units)

# Pull out the 500 hPa Heights
hght = data.variables['Geopotential_height'][:].squeeze() * units.meter
uwnd = data.variables['u_wind'][:].squeeze() * units('m/s')
vwnd = data.variables['v_wind'][:].squeeze() * units('m/s')

# Calculate the magnitude of the wind speed in kts
sped = get_wind_speed(uwnd, vwnd).to('knots')


##################################
# Set up the projection for LCC
plotcrs = ccrs.LambertConformal(central_longitude=-100.0, central_latitude=45.0)
datacrs = ccrs.PlateCarree(central_longitude=0.)

states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lakes',
        scale='50m',
        facecolor='none')


##################################
# Subset and smooth

# Subset the data arrays to grab only 500 hPa
hght_500 = hght[lev_500]
uwnd_500 = uwnd[lev_500]
vwnd_500 = vwnd[lev_500]

# Smooth the 500-hPa geopotential height field
# Be sure to only smooth the 2D field
Z_500 = ndimage.gaussian_filter(hght_500, sigma=5, order=0)


##################################
# Plot the contours

# Start plot with new figure and axis
fig = plt.figure(figsize=(17., 11.))
ax = plt.subplot(1, 1, 1, projection=plotcrs)

# Add some titles to make the plot readable by someone else
plt.title('500-hPa Geo Heights (m; black), Smoothed 500-hPa Geo. Heights (m; red)',
          loc='left')
plt.title('VALID: %s'.format(vtimes[0]), loc='right')

# Set GAREA and add map features
ax.set_extent([-125., -67., 22., 52.], ccrs.PlateCarree())
ax.coastlines('50m', edgecolor='black', linewidth=0.75)
ax.add_feature(states_provinces, edgecolor='black', linewidth=0.5)

# Set the CINT
clev500 = np.arange(5100, 6000, 60)

# Plot smoothed 500-hPa contours
cs2 = ax.contour(lons, lats, Z_500, clev500, colors='red',
                 linewidths=3, linestyles='solid', transform=datacrs)
c2 = plt.clabel(cs2, fontsize=12, colors='red', inline=1, inline_spacing=8,
                fmt='%i', rightside_up=True, use_clabeltext=True)

# Contour the 500 hPa heights with labels
cs = ax.contour(lons, lats, hght_500, clev500, colors='black',
                linewidths=2.5, linestyles='solid', alpha=0.6, transform=datacrs)
cl = plt.clabel(cs, fontsize=12, colors='k', inline=1, inline_spacing=8,
                fmt='%i', rightside_up=True, use_clabeltext=True)

plt.show()
