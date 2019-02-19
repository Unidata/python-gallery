"""
500 hPa Geopotential Heights and Winds
======================================

Classic 500-hPa plot using NAM analysis file.

This example uses example data from the NAM anlysis for 12 UTC 31
October 2016 and uses xarray as the main read source with using Cartopy
for plotting a CONUS view of the 500-hPa geopotential heights, wind
speed, and wind barbs.

"""


######################################################################
# Import the needed modules.
#

from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
from scipy.ndimage import gaussian_filter
import xarray as xr


######################################################################
# The following code reads the example data using the xarray open_dataset
# function and prints the coordinate values that are associated with the
# various variables contained within the file.
#

ds = xr.open_dataset('https://thredds.ucar.edu/thredds/dodsC/casestudies/'
                     'python-gallery/NAM_20161031_1200.nc')
ds.coords


######################################################################
# Data Retrieval
# --------------
#
# This code retrieves the necessary data from the file and completes some
# smoothing of the geopotential height and wind fields using the SciPy
# function gaussian_filter. A nicely formated valid time (vtime) variable
# is also created.
#

# Grab lat/lon values (NAM will be 2D)
lats = ds.lat.data
lons = ds.lon.data

# Select and grab data
hght = ds['Geopotential_height_isobaric']
uwnd = ds['u-component_of_wind_isobaric']
vwnd = ds['v-component_of_wind_isobaric']

# Select and grab 500-hPa geopotential heights and wind components, smooth with gaussian_filter
hght_500 = gaussian_filter(hght.sel(isobaric=500).data[0], sigma=3.0)
uwnd_500 = gaussian_filter(uwnd.sel(isobaric=500).data[0], sigma=3.0) * units('m/s')
vwnd_500 = gaussian_filter(vwnd.sel(isobaric=500).data[0], sigma=3.0) * units('m/s')

# Use MetPy to calculate the wind speed for colorfill plot, change units to knots from m/s
sped_500 = mpcalc.wind_speed(uwnd_500, vwnd_500).to('kt')

# Create a clean datetime object for plotting based on time of Geopotential heights
vtime = datetime.strptime(str(ds.time.data[0].astype('datetime64[ms]')),
                          '%Y-%m-%dT%H:%M:%S.%f')


######################################################################
# Map Creation
# ------------
#
# This next set of code creates the plot and draws contours on a Lambert
# Conformal map centered on -100 E longitude. The main view is over the
# CONUS with geopotential heights contoured every 60 m and wind speed in
# knots every 20 knots starting at 30 kt.
#

# Set up the projection that will be used for plotting
mapcrs = ccrs.LambertConformal(central_longitude=-100,
                               central_latitude=35,
                               standard_parallels=(30, 60))

# Set up the projection of the data; if lat/lon then PlateCarree is what you want
datacrs = ccrs.PlateCarree()

# Start the figure and create plot axes with proper projection
fig = plt.figure(1, figsize=(14, 12))
ax = plt.subplot(111, projection=mapcrs)
ax.set_extent([-130, -72, 20, 55], ccrs.PlateCarree())

# Add geopolitical boundaries for map reference
ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))

# Plot 500-hPa Colorfill Wind Speeds in knots
clevs_500_sped = np.arange(30, 150, 20)
cf = ax.contourf(lons, lats, sped_500, clevs_500_sped, cmap=plt.cm.BuPu,
                 transform=datacrs)
plt.colorbar(cf, orientation='horizontal', pad=0, aspect=50)

# Plot 500-hPa Geopotential Heights in meters
clevs_500_hght = np.arange(0, 8000, 60)
cs = ax.contour(lons, lats, hght_500, clevs_500_hght, colors='black',
                transform=datacrs)
plt.clabel(cs, fmt='%d')

# Plot 500-hPa wind barbs in knots, regrid to reduce number of barbs
ax.barbs(lons, lats, uwnd_500.to('kt').m, vwnd_500.to('kt').m, pivot='middle',
         color='black', regrid_shape=20, transform=datacrs)

# Make some nice titles for the plot (one right, one left)
plt.title('500-hPa NAM Geopotential Heights (m), Wind Speed (kt),'
          ' and Wind Barbs (kt)', loc='left')
plt.title('Valid Time: {}'.format(vtime), loc='right')

# Adjust image and show
plt.subplots_adjust(bottom=0, top=1)
plt.show()
