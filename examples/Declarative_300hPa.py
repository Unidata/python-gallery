"""
MetPy Declarative - 300 hPa
===========================

By: Kevin Goebbert

This example uses the declarative syntax available through the MetPy
package to allow a more convenient method for creating simple maps of
atmospheric data. The key thing the declarative language does is to
reduce the number of packages that users will need to know in detail and
instead allow them to set key parameters to get the map they desire. One
key element is the use of xarray as the data object, which allows
coordinate information to be associated with atmospheric variables.

"""

from datetime import datetime

import metpy.calc as mpcalc
import metpy.plots as mpplots
from metpy.units import units
import xarray as xr


######################################################################
# Open dataset using xarray module and subset global GFS to be over the
# CONUS.
#

ds = xr.open_dataset('https://thredds.ucar.edu/thredds/dodsC/casestudies/python-gallery/'
                     'GFS_20101026_1200.nc').sel(lon=slice(360-150, 360-50, 2),
                                                 lat=slice(65, 20, 2))


######################################################################
# Calculate Variable and Add to Dataset
# -------------------------------------
#
# Here it is demonstrated how you can calculate a new variable and add it
# to the xarray dataset (ds) so that it can be plotted with the
# declarative syntax. The key to adding a variable to an xarray dataset
# for use in the declarative syntax is the need to add a ``grid_mapping``
# and ``units`` attribute.
#

# Calculate New Variables and place into Xarray Dataset
uwnd = ds['u-component_of_wind_isobaric']
vwnd = ds['v-component_of_wind_isobaric']

# Compute wind speed using MetPy
wspd = mpcalc.wind_speed(uwnd, vwnd)

# Place wind speed (wspd) into xarray dataset and attach needed attributes
ds = ds.assign(wind_speed=(tuple(uwnd.dims)[:4], wspd.m,
                           {'grid_mapping': uwnd.attrs['grid_mapping'],
                            'units': str(wspd.units)}))

# Convert wind speed units to knots within the dataset
ds.wind_speed.metpy.convert_units('kt')
ds['u-component_of_wind_isobaric'].metpy.convert_units('kt')
ds['v-component_of_wind_isobaric'].metpy.convert_units('kt')


######################################################################
# Declarative Plot
# ----------------
#
# The following settings create a single panel map plot of 300 hPa
# geopotential heights, wind speed, and wind barbs.
#

# Countour Plot of Geopotential Heights
contour = mpplots.ContourPlot()
contour.data = ds
contour.time = datetime(2010, 10, 31, 12)
contour.field = 'Geopotential_height_isobaric'
contour.level = 300 * units.hPa
contour.linecolor = 'black'
contour.linestyle = '-'
contour.linewidth = 2
contour.clabels = True
contour.contours = list(range(0, 20000, 120))

# Colorfilled Plot of Wind Speed
cfill = mpplots.FilledContourPlot()
cfill.data = ds
cfill.field = 'wind_speed'
cfill.level = 300 * units.hPa
cfill.colormap = 'BuPu'
cfill.contours = list(range(50, 171, 20))
cfill.colorbar = 'vertical'

# Plot wind barbs
barb = mpplots.BarbPlot()
barb.data = ds
barb.level = 300 * units.hPa
barb.field = ['u-component_of_wind_isobaric', 'v-component_of_wind_isobaric']
barb.skip = (3, 3)
barb.color = 'black'
barb.barblength = 6.5
barb.earth_relative = False

# Panel for plot with Map features
panel = mpplots.MapPanel()
panel.layout = (1, 1, 1)
panel.area = (-124, -72, 20, 53)
panel.projection = 'lcc'
panel.layers = ['coastline', 'borders', 'states', 'land']
panel.plots = [cfill, contour, barb]

# Bringing it all together
pc = mpplots.PanelContainer()
pc.size = (15, 9)
pc.panels = [panel]

pc.show()
