"""
=================
Surface Meteogram
=================

Plot a surface meteogram from METAR data available on the Unidata Thredds Server.

Plot a four panel meteogram and use MetPy calculations to calculate or modify
some variables.
"""
######################################
# Import appropriate libraries
from datetime import datetime, timedelta

from matplotlib.dates import AutoDateLocator, DateFormatter
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS

########################################
# Begin Data Ingest
# -----------------

# Request METAR data from TDS
metar = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/nws/metar/ncdecoded/catalog.xml')
dataset = list(metar.datasets.values())[0]
print(list(dataset.access_urls))

########################################
# What variables are available in dataset?

# Access netcdf subset and use siphon to request data
ncss_url = dataset.access_urls['NetcdfSubset']
ncss = NCSS(ncss_url)
print(ncss.variables)

########################################
# Set query to get desired data from Thredds server

# get current date and time
now = datetime.utcnow()
now = datetime(now.year, now.month, now.day, now.hour)

# define time range you want the data for
start = now - timedelta(days=1)
end = now

# build the query
query = ncss.query()
query.lonlat_point(-90.08, 32.32)
query.time_range(start, end)
query.variables('air_temperature', 'dew_point_temperature', 'wind_speed',
                'precipitation_amount_hourly', 'inches_ALTIM',
                'air_pressure_at_sea_level', 'wind_from_direction')
query.accept('netcdf')

# Get the netcdf dataset
data = ncss.get_data(query)
print(list(data.variables))

########################################
# Begin parsing out information from file including data

# Get the station ID
station_id = data['station_id'][:].tostring()
station_id = station_id.decode('utf-8')
print(station_id)

# Get time into a datetime object
time = [datetime.fromtimestamp(t) for t in data['time']]
time = sorted(time)
print(time)

temp = data.variables['air_temperature'][:] * units('degC')
dewp = data.variables['dew_point_temperature'][:] * units('degC')
slp = data.variables['inches_ALTIM'][:] * units('inHg')
wspd = data.variables['wind_speed'][:] * units('m/s')
wdir = data.variables['wind_from_direction'][:] * units('degree')

########################################
# Use MetPy Calculations to calculate RH
# --------------------------------------

# Get ambient partial pressure, use to calculate mixing ratio
es = mpcalc.saturation_vapor_pressure(dewp)
mixr = mpcalc.mixing_ratio(es, slp)

# Calculate vapor pressure
vp = mpcalc.vapor_pressure(slp, mixr)

# Calculate saturation vapor pressure
svp = mpcalc.saturation_vapor_pressure(temp)

# Calculate relative humidity as a percentage
rh = (vp / svp) * 100


########################################
# Make Meteogram Plot
# -------------------

# Create the plots
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(12, 10))
ax1.plot(time, temp, ls='solid', marker='o', color='tab:red', ms=5)
ax1.plot(time, dewp, ls='solid', color='tab:green', marker='*', ms=8)
ax1.fill_between(time, dewp.min().m-5, dewp.m, color='tab:green', alpha=0.5)
ax1.fill_between(time, dewp.m, temp.m, color='tab:red', alpha=0.5)
ax1.set_ylabel(r'$Temp\/(^oC)$''\n'r'$Dew\/Point\/Temp\/(^oC)$', fontsize=12)
ax1.set_ylim(dewp.min().m - 5, temp.max().m + 5)
ax1.grid(True)

ax2.bar(time, wspd, width=.01, align='center', color='skyblue')
ax2.set_ylabel(r'$Wind\/Speed\/(m/s)$', fontsize=12)
ax2b = ax2.twinx()
ax2b.plot(time, wdir, marker='d', ls='None')
ax2b.set_ylim(-20, 360)
ax2b.set_ylabel(r'$Wind\ Dir\/(^{o})$')
ax2.grid(True)

ax3.plot(time, rh, color='darkgreen', marker='^')
ax3.set_ylabel(r'$Relative\/Humidity\/(\%)$', fontsize=12)
ax3.fill_between(time, 0, rh, color='palegreen')
ax3.set_ylim(0, 100)
ax3.grid(True)

ax4.plot(time, slp, ls='--', color='brown', lw=3)
ax4.set_ylabel(r'$Pressure\/(in-Hg)$')
ax4.set_ylim(slp.min().m-0.05, slp.max().m+0.05)
ax4.grid(True)

locator = AutoDateLocator()
fmt = DateFormatter('%H:%M')

ax1.xaxis.set_major_locator(locator)
ax1.xaxis.set_major_formatter(fmt)
ax1.autoscale_view()
ax1.set_title('Site: {}     Date: {:%Y/%m/%d}'.format(station_id, time[0]), fontsize=16)
ax4.set_xlabel(r'$Hour\/of\/day$', fontsize=14)
fig.autofmt_xdate()
plt.show()
