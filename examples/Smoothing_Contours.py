
# coding: utf-8

# # Smoothing Contours
# By: Kevin Goebbert
# 
# Date: 13 April 2017

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from netCDF4 import num2date
import scipy.ndimage as ndimage
from scipy.interpolate import interp2d
from datetime import datetime
from siphon.ncss import NCSS
from metpy.units import units
from metpy.calc import get_wind_speed


# In[2]:

year = '2016'
month = '04'
day = '16'
hour = '18'
ym = year+month
ymd = year+month+day

# Set up netCDF Subset Service link
ncss = NCSS('http://nomads.ncdc.noaa.gov/thredds/ncss/grid/namanl/'+ym+'/'+ymd+'/namanl_218_'+ymd+'_'+hour+'00_000.grb')

# Data Query
hgt = ncss.query().time(datetime(int(year),int(month),int(day),int(hour)))
hgt.variables('Geopotential_height','u_wind','v_wind').add_lonlat()

# Actually getting the data
data = ncss.get_data(hgt)


# In[3]:

# Get dimension names to pull appropriate variables
dtime = data.variables['Geopotential_height'].dimensions[0]
dlev  = data.variables['Geopotential_height'].dimensions[1]
dlat  = data.variables['Geopotential_height'].dimensions[2]
dlon  = data.variables['Geopotential_height'].dimensions[3]

# Get lat and lon data, as well as time data and metadata
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
lons[lons>180] = lons[lons>180]-360

# Need 2D lat/lons for plotting, do so if necessary
if (lats.ndim < 2):
    lons, lats = np.meshgrid(lons,lats)

# Determine the level of 500 hPa
levs = data.variables[dlev][:]
lev_500 = np.where(levs==500)[0][0]

times = data.variables[dtime]
# Create more useable times for output
vtimes = num2date(times[:],times.units)

# Pull out the 500 hPa Heights
hght = data.variables['Geopotential_height'][:].squeeze() * units.meter
uwnd = data.variables['u_wind'][:].squeeze() * units('m/s')
vwnd = data.variables['v_wind'][:].squeeze() * units('m/s')

# Calculate the magnitude of the wind speed in kts
sped = get_wind_speed(uwnd, vwnd).to('knots')


# In[4]:

# Set up the projection for LCC
plotcrs = ccrs.LambertConformal(central_longitude=-100.0, central_latitude=45.0)
datacrs = ccrs.PlateCarree(central_longitude=0.)

states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lakes',
        scale='50m',
        facecolor='none')


# In[5]:

# Subset the data arrays to grab only 500 hPa
hght_500 = hght[lev_500]
uwnd_500 = uwnd[lev_500]
vwnd_500 = vwnd[lev_500]

# Smooth the 500-hPa geopotential height field
# Be sure to only smooth the 2D field
Z_500 = ndimage.gaussian_filter(hght_500, sigma=5, order=0)


# In[12]:

get_ipython().magic('matplotlib inline')
# Start plot with new figure and axis
fig = plt.figure(1,figsize=(17.,11.))
ax  = plt.subplot(111,projection=plotcrs)

# Add some titles to make the plot readable by someone else
plt.title('500-hPa Geo Heights (m; black), Smoothed 500-hPa Geo. Heights (m; red)',loc='left')
plt.title('VALID: %s' %(vtimes[0]),loc='right')

# Set GAREA and add map features
ax.set_extent([-125.,-67.,22.,52.], ccrs.PlateCarree())
ax.coastlines('50m',edgecolor='black',linewidth=0.75)
ax.add_feature(states_provinces,edgecolor='black',linewidth=0.5)


# Set the CINT
clev500 = np.arange(5100,6000,60)
# Plot smoothed 500-hPa contours
cs2 = ax.contour(lons, lats, Z_500, clev500, colors='red',
                 linewidths=3, linestyles='solid', transform=datacrs)
c2 = plt.clabel(cs2, fontsize=12, colors='red',inline=1, inline_spacing=8, 
                fmt='%i', rightside_up=True, use_clabeltext=True)

# Contour the 500 hPa heights with labels
cs = ax.contour(lons, lats, hght_500, clev500, colors='black',
                linewidths=2.5, linestyles='solid', alpha = 0.6, transform=datacrs)
cl = plt.clabel(cs, fontsize=12, colors='k',inline=1, inline_spacing=8, 
                fmt='%i', rightside_up=True, use_clabeltext=True)



# Transform Vectors before plotting, then plot wind barbs.
#wslice = slice(None, None, 15)
#ax.barbs(lons[wslice,wslice], lats[wslice,wslice],
#         uwnd_500[wslice,wslice].m, vwnd_500[wslice,wslice].m,
#         length=7, transform=datacrs)

# Save and show
#plt.savefig('500hPa_heights_'+vtimes[0].strftime('%Y%m%d_%H')+'.png',dpi=150)
plt.show()


# In[ ]:



