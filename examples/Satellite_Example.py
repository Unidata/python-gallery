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
from datetime import datetime
from urllib.request import urlopen

import cartopy.crs as ccrs
import cartopy.feature as cfeat
from matplotlib import patheffects
import matplotlib.pyplot as plt
from metpy.io import GiniFile
from metpy.plots.ctables import registry
from metpy.units import units
from netCDF4 import num2date
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS


##############################################
# Get satellite data and set projection based on that data.

# Scan the catalog and download the data
satcat = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/satellite/'
                    'WV/EAST-CONUS_4km/current/catalog.xml')
dataset = satcat.datasets[list(satcat.datasets)[0]]
gini_ds = GiniFile(urlopen(dataset.access_urls['HTTPServer'])).to_dataset()

# Pull parts out of the data file
data_var = gini_ds.variables['WV']
x = gini_ds.variables['x'][:]
y = gini_ds.variables['y'][:]
time_var = gini_ds.variables['time']
proj_var = gini_ds.variables[data_var.grid_mapping]

# Set up the projection based on satellite projection
globe = ccrs.Globe(ellipse='sphere', semimajor_axis=proj_var.earth_radius,
                   semiminor_axis=proj_var.earth_radius)

proj = ccrs.LambertConformal(central_longitude=proj_var.longitude_of_central_meridian,
                             central_latitude=proj_var.latitude_of_projection_origin,
                             standard_parallels=[proj_var.standard_parallel],
                             globe=globe)

##############################################
# Use Siphon to obtain data that is close to the time of the satellite file

gfscat = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/'
                    'NCEP/GFS/Global_0p5deg/catalog.xml')
dataset = gfscat.datasets[list(gfscat.datasets)[1]]
ncss = NCSS(dataset.access_urls['NetcdfSubset'])
now = datetime.utcnow()

query = ncss.query()
query.variables('Geopotential_height_isobaric',
                'u-component_of_wind_isobaric',
                'v-component_of_wind_isobaric')
query.add_lonlat().vertical_level(300*100)
query.time(now)
query.lonlat_box(north=65, south=15, east=310, west=220)
data = ncss.get_data(query)

##############################################
# Pull out specific variables and attach units.

hght_300 = data.variables['Geopotential_height_isobaric'][:].squeeze() * units.meter
uwnd_300 = data.variables['u-component_of_wind_isobaric'][:].squeeze() * units('m/s')
vwnd_300 = data.variables['v-component_of_wind_isobaric'][:].squeeze() * units('m/s')

Z_300 = ndimage.gaussian_filter(hght_300, sigma=4, order=0)

lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
time = data.variables[data.variables['Geopotential_height_isobaric'].dimensions[0]]
vtime = num2date(time[:], time.units)

##############################################
# Create figure with an overlay of WV Imagery with 300-hPa Heights and Wind

# Create the figure
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=proj)

# Add mapping information
ax.coastlines(resolution='50m', color='black')
state_boundaries = cfeat.NaturalEarthFeature(category='cultural',
                                             name='admin_1_states_provinces_lakes',
                                             scale='50m', facecolor='none')

ax.add_feature(state_boundaries, edgecolor='black', linestyle=':')
ax.add_feature(cfeat.BORDERS, linewidth='2', edgecolor='black')

# Plot the image with our colormapping choices
wv_norm, wv_cmap = registry.get_with_steps('WVCIMSS', 0, 1)
im = ax.imshow(data_var[:], extent=(x[0], x[-1], y[0], y[-1]), origin='upper',
               cmap=wv_cmap, norm=wv_norm)

# Add the text, complete with outline
timestamp = num2date(time_var[:].squeeze(), time_var.units)
text = ax.text(0.99, 0.01, timestamp.strftime('%d %B %Y %H%MZ'),
               horizontalalignment='right', transform=ax.transAxes,
               color='white', fontsize='x-large', weight='bold')
text.set_path_effects([patheffects.withStroke(linewidth=2, foreground='black')])

# PLOT 300-hPa Geopotential Heights and Wind Barbs
ax.set_extent([-112, -65, 20, 59], ccrs.Geodetic())
cs = ax.contour(lon, lat, Z_300, colors='black', transform=ccrs.PlateCarree())
ax.clabel(cs, fontsize=12, colors='k', inline=1, inline_spacing=8,
          fmt='%i', rightside_up=True, use_clabeltext=True)

ax.barbs(lon, lat, uwnd_300.to('knots').m, vwnd_300.to('knots').m, color='tab:red',
         length=7, regrid_shape=15, pivot='middle', transform=ccrs.PlateCarree())

ax.set_title('300-hPa Geo Heights (m; black) and Wind Barbs (kt)', loc='left')
ax.set_title('Valid: {}'.format(vtime[0]), loc='right')

plt.show()
