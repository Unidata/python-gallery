"""
======================
GOES ABI + GLM example
======================

This example plots the most recent GOES-16 data with an overlay of the current
GOES-16 Geostationary Lightning Mapper flash extent density.
"""

from datetime import datetime
import json
import urllib.request

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import patheffects
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import metpy  # noqa: F401
import metpy.calc as mpcalc
from metpy.plots.ctables import registry
from metpy.units import units
import numpy as np
from siphon.catalog import TDSCatalog
import xarray as xr
from xarray.backends import NetCDF4DataStore

# e.g. OR_GLM-L2-GLMC-M3_G16_s20191930223000_e20191930224000_c20191930224350.nc
# leave off last char (13 not 14), to leave off 1/10 s part
goes_start_regex =  r'_s(?P<strptime>\d{13})'
goes_strptime = '%Y%j%H%M%S'

def get_glm_image(date=datetime.utcnow(), region='CONUS'):
    """Return dataset of GOES-16 GLM data."""
    cat = TDSCatalog('https://thredds-test.unidata.ucar.edu/thredds/catalog/satellite/goes/'
                     'east/products/GeostationaryLightningMapper/{}/{:%Y%m%d}/'
                     'catalog.xml'.format(region, date))
    ds = cat.datasets.filter_time_nearest(date, regex=goes_start_regex,
                strptime=goes_strptime)
    ds = ds.remote_access(service='CdmRemote', use_xarray=True)
    return ds


def get_goes_image(date=datetime.utcnow(), channel=8, region='CONUS'):
    """Return dataset of GOES-16 data."""
    cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/satellite/goes/east/products/'
                     'CloudAndMoistureImagery/{}/Channel{:02d}/{:%Y%m%d}/'
                     'catalog.xml'.format(region, channel, date))
    ds = cat.datasets.filter_time_nearest(date, regex=goes_start_regex,
                strptime=goes_strptime)
    ds = ds.remote_access(service='OPENDAP')
    # ••• Change this to xarray direct, too? •••
    ds = NetCDF4DataStore(ds)
    ds = xr.open_dataset(ds)
    return ds

# lat_lon_box = [-124.5, -105, 38.5, 50]

ds = get_goes_image(channel=14)
ds_glm = get_glm_image()

# Parse out the projection data from the satellite files
dat = ds.metpy.parse_cf('Sectorized_CMI')
proj = dat.metpy.cartopy_crs
dat_glm = ds_glm.metpy.parse_cf('flash_extent_density')
proj_glm = dat_glm.metpy.cartopy_crs

fig = plt.figure(figsize=(1.375 * 40, 40))
ax = fig.add_subplot(1, 1, 1, projection=proj)

def plot_image(dat, ax, cmap='WVCIMSS_r', cmap_range=(195,265)):
    # Pull out what we need from the GOES netCDF file
    x = dat['x']
    y = dat['y']

    # Make the plot
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

    try:
        mpl_norm, mpl_cmap = registry.get_with_range(cmap, *cmap_range)
    except KeyError:
        # If it's not in the MetPy registry, revert to ordinary matplotlib
        mpl_norm, mpl_cmap = Normalize(*cmap_range), cmap 

    im = ax.imshow(dat, extent=(x.min(), x.max(), y.min(), y.max()),
                   origin='upper')
    im.set_cmap(mpl_cmap)
    im.set_norm(mpl_norm)

plot_image(dat, ax, cmap='gray', cmap_range=(-70+273.5,40+273.5))
plot_image(dat_glm, ax, cmap='YlOrRd', cmap_range=(0,32))           

ax.add_feature(cfeature.BORDERS, linewidth=8, edgecolor='black')
ax.add_feature(cfeature.COASTLINE, linewidth=8, edgecolor='black')
ax.add_feature(cfeature.STATES.with_scale('50m'), linestyle='-',
               edgecolor='black', linewidth=4)

timestamp = datetime.strptime(ds.start_date_time, '%Y%j%H%M%S')
timestamp_glm = datetime.strptime(ds_glm.time_coverage_start, '%Y-%m-%dT%H:%M:%SZ')
abi_label = 'ABI Ch 14: '+timestamp.strftime('%d %B %Y %H%MZ')
glm_label = 'GLM FED: '+timestamp_glm.strftime('%d %B %Y %H%MZ')

text_time = ax.text(0.01, 0.01, abi_label,
                    horizontalalignment='left', transform=ax.transAxes,
                    color='white', fontsize=72, weight='bold')
text_time_glm = ax.text(0.01, 0.05, glm_label,
                    horizontalalignment='left', transform=ax.transAxes,
                    color='white', fontsize=72, weight='bold')

outline_effect = [patheffects.withStroke(linewidth=15, foreground='black')]
text_time.set_path_effects(outline_effect)
text_time_glm.set_path_effects(outline_effect)

# ax.set_extent(lat_lon_box)
ax.gridlines(linestyle=':', color='black', linewidth=2)

plt.savefig('sat_image.png')
