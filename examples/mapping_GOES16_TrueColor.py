"""
==========================
GOES-16: True Color Recipe
==========================
By: [Brian Blaylock](http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/home.html)
with help from Julien Chastang (UCAR-Unidata).

Additional notebooks analyzing GOES-16 and other data can be found in [Brian's
GitHub repository](https://github.com/blaylockbk/pyBKB_v3/).

This notebook shows how to make a true color image from the GOES-16
Advanced Baseline Imager (ABI) level 2 data. We will plot the image with
matplotlib and Cartopy. The image can be displayed on any map projection after
applying a transformation using pyproj. The methods shown here are stitched
together from the following useful information found online:


- [**True Color RGB Recipe**](http://cimss.ssec.wisc.edu/goes/OCLOFactSheetPDFs/ABIQuickGuide_CIMSSRGB_v2.pdf)
- [ABI Bands Quick Information Guides](https://www.goes-r.gov/education/ABI-bands-quick-info.html)
- [Open Commons Consortium](http://edc.occ-data.org/goes16/python/)
- [GeoNetCast Blog](https://geonetcast.wordpress.com/2017/07/25/geonetclass-manipulating-goes-16-data-with-python-part-vi/)
- [Proj documentation](https://proj4.org/operations/projections/geos.html?highlight=geostationary)
- [Pyproj documentation](http://jswhit.github.io/pyproj/pyproj.Proj-class.html)

True color images are an RGB composite of the following three channels:

|        --| Wavelength   | Channel | Description   |
|----------|--------------|---------|---------------|
| **Red**  | 0.64 &#181;m |    2    | Red Visible   |
| **Green**| 0.86 &#181;m |    3    | Veggie Near-IR|
| **Blue** | 0.47 &#181;m |    1    | Blue Visible  |

For this demo, you will need GOES-16 ABI level 2 data. You can get GOES-16 files
from NOAA's GOES archive on [Amazon
S3](https://aws.amazon.com/public-datasets/goes/). I created a [web
interface](http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi?domain=C&product=ABI-L2-MCMIP&hour=0)
to easily download files from the Amazon archive. For scripted or bulk
downloads, you should use `rclone` or `AWS CLI`. You may also download files
from the [Environmental Data Commons](http://edc.occ-data.org/goes16/getdata/)
and [NOAA
CLASS](https://www.avl.class.noaa.gov/saa/products/search?sub_id=0&datatype_family=GRABIPRD&submit.x=25&submit.y=9).

This example uses the **level 2 _multiband_ formatted file for the _CONUS_
domain**
([ABI-L2-MCMIPC](http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi?domain=C&product=ABI-L2-MCMIP&hour=0)).
The multiband file is easiest to work with becuase it contains all 16 channels
on the same 2 km grid. However, some channels have higher resolution. Plotting
the full resolution images will take some additional work, not described
here. Specifically, you will have to download three separate files, one for each
channel, and subsample the red channel 0.5 km grid to a 1 km grid.

I previously downloaded the following file from Amazon Web Services

OR_ABI-L2-MCMIPC-M3_G16_s20181781922189_e20181781924562_c20181781925075.nc

OR     - Indicates the system is operational
ABI    - Instrument type
L2     - Level 2 Data
MCMIP  - Multichannel Cloud and Moisture Imagery products
c      - CONUS file (created every 5 minutes).
M3     - Scan mode
G16    - GOES-16
s##### - Scan start: 4 digit year, 3 digit day of year (Julian day),
hour, minute, second, tenth second
e##### - Scan end
c##### - File Creation
.nc    - NetCDF file extension

"""  # noqa: E501

######################################################################
# First, import the libraries we will use

from datetime import datetime

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import metpy  # noqa: F401
import numpy as np
import xarray


######################################################################
# Open the GOES-16 NetCDF File
# ----------------------------

# Using xarray, I assign the opened file to the variable C for the CONUS domain.

FILE = ('http://ramadda-jetstream.unidata.ucar.edu/repository/opendap'
        '/4ef52e10-a7da-4405-bff4-e48f68bb6ba2/entry.das#fillmismatch')
C = xarray.open_dataset(FILE)

######################################################################
# Date and Time Information
# ----------------------------
# Each file represents the data collected during one scan sequence for the
# domain. There are several different time stamps in this file, which are also
# found in the file's name.

# Scan's start time, converted to datetime object
scan_start = datetime.strptime(C.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

# Scan's end time, converted to datetime object
scan_end = datetime.strptime(C.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')

# File creation time, convert to datetime object
file_created = datetime.strptime(C.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')

# The 't' variable is the scan's midpoint time
# I'm not a fan of numpy datetime, so I convert it to a regular datetime object
midpoint = str(C['t'].data)[:-8]
scan_mid = datetime.strptime(midpoint, '%Y-%m-%dT%H:%M:%S.%f')

print('Scan Start    : {}'.format(scan_start))
print('Scan midpoint : {}'.format(scan_mid))
print('Scan End      : {}'.format(scan_end))
print('File Created  : {}'.format(file_created))
print('Scan Duration : {:.2f} minutes'.format((scan_end-scan_start).seconds/60))


######################################################################
# True Color Recipe
# -----------------
#
# Color images are a Red-Green-Blue (RGB) composite of three different
# channels. We will assign the following channels as our RGB values:
#
# | --                   | RED         | GREEN          | BLUE         |
# |----------------------|-------------|----------------|--------------|
# | **Name**             | Red Visible | Near-IR Veggie | Blue Visible |
# | **Wavelength**       | 0.64 µm     | 0.86 µm        | 0.47 µm      |
# | **Channel**          | 2           | 3              | 1            |
# | **Units**            | Reflectance | Reflectance    | Reflectance  |
# | **Range of Values**  | 0-1         | 0-1            | 0-1          |
# | **Gamma Correction** | 2.2         | 2.2            | 2.2          |
#
# RGB values must be between 0 and 1, same as the range of values of the
# reflectance channels. A gamma correction is applied to control the brightness
# and make the image not look too dark where `corrected_value =
# value^(1/gamma)`. Most displays have a decoding gamma of 2.2
# ([source1](https://en.wikipedia.org/wiki/Gamma_correction),
# [source2](https://www.cambridgeincolour.com/tutorials/gamma-correction.htm)).
#
# The GREEN "veggie" channel on GOES-16 does not measure visible green
# light. Instead, it measures a near-infrared band sensitive to chlorophyll. We
# could use that channel in place of green, but it would make the green in our
# image appear too vibrant. Instead, we will tone-down the green channel by
# interpolating the value to simulate a natural green color.
#
# \begin{equation}
# TrueGreen = (0.48358168*RED) + (0.06038137*GREEN) + (0.45706946*BLUE)
# \end{equation}
#
# or, a simple alternative ([CIMSS Natural True
# Color](http://cimss.ssec.wisc.edu/goes/OCLOFactSheetPDFs/ABIQuickGuide_CIMSSRGB_v2.pdf)):
#
# \begin{equation}
# TrueGreen = (0.45*RED) + (0.1*GREEN) + (0.45*BLUE)
# \end{equation}
#
# The multiband formatted file we loaded is convenient becuase all the GOES
# channels are in the same NetCDF file. Next, we will assign our variables R, G,
# and B as the data for each channel.

# Confirm that each band is the wavelength we are interested in
for band in [2, 3, 1]:
    print('{} is {:.2f} {}'.format(
        C['band_wavelength_C{:02d}'.format(band)].long_name,
        float(C['band_wavelength_C{:02d}'.format(band)][0]),
        C['band_wavelength_C{:02d}'.format(band)].units))

######################################################################

# Load the three channels into appropriate R, G, and B variables
R = C['CMI_C02'].data
G = C['CMI_C03'].data
B = C['CMI_C01'].data

######################################################################

# Apply range limits for each channel. RGB values must be between 0 and 1
R = np.clip(R, 0, 1)
G = np.clip(G, 0, 1)
B = np.clip(B, 0, 1)

######################################################################

# Apply a gamma correction to the image
gamma = 2.2
R = np.power(R, 1/gamma)
G = np.power(G, 1/gamma)
B = np.power(B, 1/gamma)

######################################################################

# Calculate the "True" Green
G_true = 0.45 * R + 0.1 * G + 0.45 * B
G_true = np.maximum(G_true, 0)
G_true = np.minimum(G_true, 1)

######################################################################
# Simple Image
# -----------------
#
# Use `plt.imshow` to get a quick look at the channels and RGB composite we
# created.
#
# First, we plot each channel individually. The deeper the color means the
# satellite is observing more light in that channel. Clouds appear white becuase
# they reflect lots of red, green, and blue light. You will also notice that the
# land reflects a lot of "green" in the veggie channel becuase this channel is
# sensitive to the chlorophyll.

fig, ([ax1, ax2, ax3, ax4]) = plt.subplots(1, 4, figsize=(16, 3))

ax1.imshow(R, cmap='Reds', vmax=1, vmin=0)
ax1.set_title('Red', fontweight='semibold')
ax1.axis('off')

ax2.imshow(G, cmap='Greens', vmax=1, vmin=0)
ax2.set_title('Veggie', fontweight='semibold')
ax2.axis('off')

ax3.imshow(G_true, cmap='Greens', vmax=1, vmin=0)
ax3.set_title('"True" Green', fontweight='semibold')
ax3.axis('off')

ax4.imshow(B, cmap='Blues', vmax=1, vmin=0)
ax4.set_title('Blue', fontweight='semibold')
ax4.axis('off')

plt.subplots_adjust(wspace=.02)

######################################################################
# The addition of the three channels results in a color image. We combine the
# three channels in a stacked array and display the image with `imshow` again.

# The RGB array with the raw veggie band
RGB_veggie = np.dstack([R, G, B])

# The RGB array for the true color image
RGB = np.dstack([R, G_true, B])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# The RGB using the raw veggie band
ax1.imshow(RGB_veggie)
ax1.set_title('GOES-16 RGB Raw Veggie', fontweight='semibold', loc='left',
              fontsize=12)
ax1.set_title('{}'.format(scan_start.strftime('%d %B %Y %H:%M UTC ')),
              loc='right')
ax1.axis('off')

# The RGB for the true color image
ax2.imshow(RGB)
ax2.set_title('GOES-16 RGB True Color', fontweight='semibold', loc='left',
              fontsize=12)
ax2.set_title('{}'.format(scan_start.strftime('%d %B %Y %H:%M UTC ')),
              loc='right')
ax2.axis('off')

######################################################################
# Plot with `Cartopy`:  Geostationary Projection
# ----------------------------------------------
#
# The image above is not georeferenced. You can see the land and oceans, but we
# do have enough information to draw state and country boundaries. We will use
# the `metpy.io` package to obtain the projection information. We will then use
# Cartopy to plot the image. The image is in a [geostationary projection
# ](https://proj4.org/operations/projections/geos.html?highlight=geostationary).

# We'll use the `CMI_C02` variable as a 'hook' to get the CF metadata.
dat = C.metpy.parse_cf('CMI_C02')

geos = dat.metpy.cartopy_crs

x = dat.x
y = dat.y

######################################################################
# The geostationary projection is the easiest way to plot the image on a
# map. Essentially, we are stretching the image across a map with the same
# projection and dimensions as the data.

fig = plt.figure(figsize=(15, 12))

ax = fig.add_subplot(1, 1, 1, projection=geos)

ax.imshow(np.flipud(RGB), origin='lower',
          extent=(x.min(), x.max(), y.min(), y.max()), transform=geos)

ax.coastlines(resolution='50m', color='black', linewidth=0.25)
ax.add_feature(ccrs.cartopy.feature.STATES, linewidth=0.25)

plt.title('GOES-16 True Color', loc='left', fontweight='semibold', fontsize=15)
plt.title('%s'.format(scan_start.strftime('%d %B %Y %H:%M UTC ')), loc='right')

plt.show()

######################################################################
# Using other projections
# ----------------------------------------------
#
# Changing the projections with cartopy is straightforward. Here we use
# the Lambert Conformal projection to display the GOES-16 data.

fig = plt.figure(figsize=(15, 12))

lc = ccrs.LambertConformal(central_longitude=-97.5, standard_parallels=(38.5,
                                                                        38.5))

ax = fig.add_subplot(1, 1, 1, projection=lc)
ax.set_extent([-135, -60, 10, 65], crs=ccrs.PlateCarree())

ax.imshow(np.flipud(RGB), origin='lower',
          extent=(x.min(), x.max(), y.min(), y.max()),
          transform=geos,
          interpolation='none')
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(ccrs.cartopy.feature.STATES, linewidth=0.5)

plt.title('GOES-16 True Color', loc='left', fontweight='semibold', fontsize=15)
plt.title('{}'.format(scan_start.strftime('%d %B %Y %H:%M UTC ')), loc='right')

plt.show()

######################################################################
# Plot with `Cartopy`: Plate Carrée  Cylindrical Projection
# ---------------------------------------------------------
#
# It is often useful to zoom on a specific location. This image will zoom in on
# Utah.

fig = plt.figure(figsize=(8, 8))

pc = ccrs.PlateCarree()

ax = fig.add_subplot(1, 1, 1, projection=pc)
ax.set_extent([-114.75, -108.25, 36, 43], crs=pc)

ax.imshow(np.flipud(RGB), origin='lower',
          extent=(x.min(), x.max(), y.min(), y.max()),
          transform=geos,
          interpolation='none')
ax.coastlines(resolution='50m', color='black', linewidth=1)
ax.add_feature(ccrs.cartopy.feature.STATES)

plt.title('GOES-16 True Color', loc='left', fontweight='semibold', fontsize=15)
plt.title('{}'.format(scan_start.strftime('%d %B %Y %H:%M UTC ')), loc='right')

plt.show()

######################################################################
# Overlay Nighttime IR when dark
# ------------------------------
#
# At nighttime, the visible wavelengths cannot measure anything. The entire
# domain is black. But there is information from other channels we can use to
# see where the clouds are at night. To view clouds in portions of the domain
# experiencing nighttime, we will overlay the clean infrared (IR) channel over
# the true color image.
#
# First, open a file where the scan shows partial night area and create the true
# color RGB as before.

# A GOES-16 file with half day and half night

FILE = ('http://ramadda-jetstream.unidata.ucar.edu/repository/opendap'
        '/85da3304-b910-472b-aedf-a6d8c1148131/entry.das#fillmismatch')
C = xarray.open_dataset(FILE)

# Scan's start time, converted to datetime object
scan_start = datetime.strptime(C.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

# Load the three channels into appropriate R, G, and B
R = C['CMI_C02'].data
G = C['CMI_C03'].data
B = C['CMI_C01'].data

# Apply range limits for each channel. RGB values must be between 0 and 1
R = np.clip(R, 0, 1)
G = np.clip(G, 0, 1)
B = np.clip(B, 0, 1)

# Apply the gamma correction
gamma = 2.2
R = np.power(R, 1/gamma)
G = np.power(G, 1/gamma)
B = np.power(B, 1/gamma)

# Calculate the "True" Green
G_true = 0.45 * R + 0.1 * G + 0.45 * B
G_true = np.maximum(G_true, 0)
G_true = np.minimum(G_true, 1)

# The final RGB array :)
RGB = np.dstack([R, G_true, B])

######################################################################
# _**Load the Clear IR  10.3 &#181;m channel (Band 13)**_
# -------------------------------------------------------
#
# Clean IR has units of Kelvin, so we need to normalize the temperature array
# between a range of values.

C['CMI_C13']

######################################################################
#
# Notice that the unit of the clean IR channel is *brightness temperature*, NOT
# reflectance. We will apply some bounds and a normalization becuase we must
# values between 0 and 1.

cleanIR = C['CMI_C13'].data

# Normalize the channel between a range. e.g. cleanIR =
# (cleanIR-minimum)/(maximum-minimum)
cleanIR = (cleanIR-90)/(313-90)

# Apply range limits for each channel. RGB values must be between 0 and 1
cleanIR = np.clip(cleanIR, 0, 1)

# Invert colors so that cold clouds are white
cleanIR = 1 - cleanIR

# Lessen the brightness of the coldest clouds so they don't appear so bright
# when we overlay it on the true color image
cleanIR = cleanIR/1.4

######################################################################
# Show the true color and clean IR images
# ---------------------------------------
#
# We want to overlay these two images, so the clean IR fills in the night sky on
# the True Color image. This way we can still see the clouds at night.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
ax1.imshow(np.dstack([R, G_true, B]))
ax2.imshow(np.dstack([cleanIR, cleanIR, cleanIR]))
ax1.set_title('True Color', fontweight='semibold')
ax2.set_title('Clean IR', fontweight='semibold')

######################################################################
#
# To fill in the dark area on the true color image, we will set each RGB channel
# to equal the maximum value between the visible channels and the IR
# channels. You should note that if the clean IR has really bright, cold clouds
# in the daylight, they will replace the color values in the true color image
# making the clouds appear whiter. Still, it makes a nice plot and allows you to
# see clouds when it is nigh time.

RGB_IR = np.dstack([np.maximum(R, cleanIR), np.maximum(G_true, cleanIR),
                    np.maximum(B, cleanIR)])

######################################################################

fig = plt.figure(figsize=(15, 12))

ax = fig.add_subplot(1, 1, 1, projection=geos)

ax.imshow(np.flipud(RGB_IR), origin='lower',
          extent=(x.min(), x.max(), y.min(), y.max()),
          transform=geos)

ax.coastlines(resolution='50m', color='black', linewidth=2)
ax.add_feature(ccrs.cartopy.feature.STATES)

plt.title('GOES-16 True Color and Night IR', loc='left', fontweight='semibold',
          fontsize=15)
plt.title('{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y'), loc='right'))

plt.show()

######################################################################
# Adjust Image Contrast
# ---------------------
#
# I think the color looks a little dull. We can make the colors pop out by
# adjusting the contrast. Adjusting image contrast is easy to do in Photoshop,
# and also easy to do in Python.
#
# Note: you should adjust the contrast _before_ you add in the Clean IR channel.


def contrast_correction(color, contrast):
    """
    Modify the contrast of an RGB
    See:
    https://www.dfstudios.co.uk/articles/programming/image-programming-algorithms/image-processing-algorithms-part-5-contrast-adjustment/

    Input:
        C - contrast level
    """
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.minimum(COLOR, 1)
    COLOR = np.maximum(COLOR, 0)
    return COLOR


# Amount of contrast
contrast_amount = 105

# Apply contrast correction
RGB_contrast = contrast_correction(np.dstack([R, G_true, B]), contrast_amount)

# Add in clean IR
RGB_contrast_IR = np.dstack([np.maximum(RGB_contrast[:, :, 0], cleanIR),
                             np.maximum(RGB_contrast[:, :, 1], cleanIR),
                             np.maximum(RGB_contrast[:, :, 2], cleanIR)])

######################################################################

fig = plt.figure(figsize=(15, 12))

ax1 = fig.add_subplot(1, 2, 1, projection=geos)
ax2 = fig.add_subplot(1, 2, 2, projection=geos)

plt.sca(ax1)
ax1.imshow(np.flipud(RGB_IR), origin='lower',
           extent=(x.min(), x.max(), y.min(), y.max()),
           transform=geos)
ax1.coastlines(resolution='50m', color='black', linewidth=2)
ax1.add_feature(ccrs.cartopy.feature.BORDERS)
plt.title('True Color and Night IR')

plt.sca(ax2)
ax2.imshow(np.flipud(RGB_contrast_IR), origin='lower',
           extent=(x.min(), x.max(), y.min(), y.max()),
           transform=geos)
ax2.coastlines(resolution='50m', color='black', linewidth=2)
ax2.add_feature(ccrs.cartopy.feature.BORDERS)
plt.title('Contrast={}'.format(contrast_amount))

plt.subplots_adjust(wspace=.02)

######################################################################
# Can we make plots for a Mesoscale scan?
# ---------------------------------------
#
# Yes. Yes we can.

# M1 is for the Mesoscale1 NetCDF file

FILE = ('http://ramadda-jetstream.unidata.ucar.edu/repository/opendap'
        '/5e02eafa-5cee-4d00-9f58-6e201e69b014/entry.das#fillmismatch')
M1 = xarray.open_dataset(FILE)

# Load the RGB arrays
R = M1['CMI_C02'][:].data
G = M1['CMI_C03'][:].data
B = M1['CMI_C01'][:].data

# Apply range limits for each channel. RGB values must be between 0 and 1
R = np.clip(R, 0, 1)
G = np.clip(G, 0, 1)
B = np.clip(B, 0, 1)

# Apply the gamma correction
gamma = 2.2
R = np.power(R, 1/gamma)
G = np.power(G, 1/gamma)
B = np.power(B, 1/gamma)

# Calculate the "True" Green
G_true = 0.45 * R + 0.1 * G + 0.45 * B
G_true = np.maximum(G_true, 0)
G_true = np.minimum(G_true, 1)

# The final RGB array :)
RGB = np.dstack([R, G_true, B])

# Scan's start time, converted to datetime object
scan_start = datetime.strptime(M1.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

# We'll use the `CMI_C02` variable as a 'hook' to get the CF metadata.
dat = M1.metpy.parse_cf('CMI_C02')

x = dat.x
y = dat.y

######################################################################

fig = plt.figure(figsize=(10, 8))

ax = fig.add_subplot(1, 1, 1, projection=lc)
ax.set_extent([-125, -70, 25, 50], crs=ccrs.PlateCarree())

ax.imshow(np.flipud(RGB), origin='lower',
          extent=(x.min(), x.max(), y.min(), y.max()),
          transform=geos)

ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(ccrs.cartopy.feature.STATES, linewidth=0.5)
ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5)

plt.title('GOES-16 True Color', fontweight='semibold', fontsize=15)
plt.title('{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y'), loc='right'))
plt.title('Mesoscale Section 1', loc='left')

plt.show()

######################################################################

fig = plt.figure(figsize=(15, 12))

ax = fig.add_subplot(1, 1, 1, projection=geos)

ax.imshow(np.flipud(RGB), origin='lower',
          extent=(x.min(), x.max(), y.min(), y.max()),
          transform=geos)

ax.coastlines(resolution='50m', color='black', linewidth=0.25)
ax.add_feature(ccrs.cartopy.feature.STATES, linewidth=0.25)

plt.title('GOES-16 True Color', fontweight='semibold', fontsize=15)
plt.title('{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y'), loc='right'))
plt.title('Mesoscale Section 1', loc='left')

plt.show()

######################################################################
# Can we do this for a Full Disk Scan? It's possible...
# -----------------------------------------------------
#
# but data files are so large that plotting is very slow. I don't need to do
# this, so I won't worry much about it. Feel free to experiment.

FILE = ('http://ramadda-jetstream.unidata.ucar.edu/repository/opendap'
        '/deb91f58-f997-41a3-a077-987529bf02b3/entry.das#fillmismatch')
F = xarray.open_dataset(FILE)

# Load the RGB arrays
R = F['CMI_C02'][:].data
G = F['CMI_C03'][:].data
B = F['CMI_C01'][:].data

# Apply range limits for each channel. RGB values must be between 0 and 1
R = np.clip(R, 0, 1)
G = np.clip(G, 0, 1)
B = np.clip(B, 0, 1)

# Apply the gamma correction
gamma = 2.2
R = np.power(R, 1/gamma)
G = np.power(G, 1/gamma)
B = np.power(B, 1/gamma)

# Calculate the "True" Green
G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G
G_true = np.maximum(G_true, 0)
G_true = np.minimum(G_true, 1)

# The final RGB array :)
RGB = np.dstack([R, G_true, B])

# We'll use the `CMI_C02` variable as a 'hook' to get the CF metadata.
dat = F.metpy.parse_cf('CMI_C02')

x = dat.x
y = dat.y

######################################################################
# _**Geostationary projection is easy**_

fig = plt.figure(figsize=(10, 8))

ax = fig.add_subplot(1, 1, 1, projection=geos)

ax.imshow(np.flipud(RGB), origin='lower',
          extent=(x.min(), x.max(), y.min(), y.max()),
          transform=geos)

ax.coastlines(resolution='50m', color='black', linewidth=1)
ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=1)

plt.title('GOES-16 True Color', fontweight='semibold', fontsize=15, loc='left')
plt.title('Full Disk\n{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
          loc='right')

plt.show()
