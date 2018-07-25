"""
============================
WV Satellite Overlay Example
============================

Plot a Gini Satellite file and overlay GFS-based data.

Using the Gini read capability of MetPy with Siphon to bring in the best GFS
data according to the current time, plot an overlay of WV imagery with 300-hPa
Geopotential Heights and Wind Barbs.
"""
##############################################
# Begin with imports, need a lot for this task.

# A whole bunch of imports
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import patheffects
import matplotlib.pyplot as plt
from metpy.io import GiniFile
from metpy.plots.ctables import registry
from metpy.units import units
from netCDF4 import num2date
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
import xarray as xr


##############################################
# Get satellite data and set projection based on that data.

# Scan the catalog and download the data
satcat = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/satellite/'
                    'WV/WEST-CONUS_4km/current/catalog.xml')
dataset = satcat.datasets[0]
f = GiniFile(dataset.remote_open())
gini_ds = xr.open_dataset(f)

# Pull parts out of the data file
dat = gini_ds.metpy.parse_cf('WV')
data_var = gini_ds.variables['WV']
x = gini_ds.variables['x'][:]
y = gini_ds.variables['y'][:]
timestamp = f.prod_desc.datetime

##############################################
# Use Siphon to obtain data that is close to the time of the satellite file

gfscat = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/'
                    'NCEP/GFS/Global_0p5deg/catalog.xml')
dataset = gfscat.datasets['Best GFS Half Degree Forecast Time Series']
ncss = dataset.subset()

query = ncss.query()
query.variables('Geopotential_height_isobaric',
                'u-component_of_wind_isobaric',
                'v-component_of_wind_isobaric')
query.add_lonlat().vertical_level(300 * 100)
query.time(timestamp)  # Use the time from the GINI file
query.lonlat_box(north=65, south=15, east=310, west=220)
data = ncss.get_data(query)

##############################################
# Pull out specific variables and attach units.

hght_300 = data.variables['Geopotential_height_isobaric'][:].squeeze() * units.meter
uwnd_300 = units('m/s') * data.variables['u-component_of_wind_isobaric'][:].squeeze()
vwnd_300 = units('m/s') * data.variables['v-component_of_wind_isobaric'][:].squeeze()

Z_300 = ndimage.gaussian_filter(hght_300, sigma=4, order=0)

lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
time = data.variables[data.variables['Geopotential_height_isobaric'].dimensions[0]]
vtime = num2date(time[:], time.units)

##############################################
# Create figure with an overlay of WV Imagery with 300-hPa Heights and Wind

# Create the figure
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=dat.metpy.cartopy_crs)

# Add mapping information
ax.coastlines(resolution='50m', color='black')
ax.add_feature(cfeature.STATES, linestyle=':')
ax.add_feature(cfeature.BORDERS, linewidth=2)

# Plot the image with our colormapping choices
wv_norm, wv_cmap = registry.get_with_range('WVCIMSS', 100, 260)
im = ax.imshow(data_var[:], extent=(x[0], x[-1], y[0], y[-1]), origin='upper',
               cmap=wv_cmap, norm=wv_norm)

# Add the text, complete with outline
text = ax.text(0.99, 0.01, timestamp.strftime('%d %B %Y %H%MZ'),
               horizontalalignment='right', transform=ax.transAxes,
               color='white', fontsize='x-large', weight='bold')
text.set_path_effects([patheffects.withStroke(linewidth=2, foreground='black')])

# PLOT 300-hPa Geopotential Heights and Wind Barbs
ax.set_extent([-132, -95, 25, 47], ccrs.Geodetic())
cs = ax.contour(lon, lat, Z_300, colors='black', transform=ccrs.PlateCarree())
ax.clabel(cs, fontsize=12, colors='k', inline=1, inline_spacing=8,
          fmt='%i', rightside_up=True, use_clabeltext=True)

ax.barbs(lon, lat, uwnd_300.to('knots').m, vwnd_300.to('knots').m, color='tab:red',
         length=7, regrid_shape=15, pivot='middle', transform=ccrs.PlateCarree())

ax.set_title('300-hPa Geo Heights (m; black) and Wind Barbs (kt)', loc='left')
ax.set_title('Valid: {}'.format(vtime[0]), loc='right')

plt.show()
