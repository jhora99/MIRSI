# # MIRSIphot - MIRSI image photometry 
# Joseph L. Hora
# Center for Astrophysics | Harvard & Smithsonian
# jhora@cfa.harvard.edu

# Import all necessary packages and functions

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import re
import os
import sys

from astropy.io import fits

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io.fits import getheader
from astropy.io.fits import getdata
from matplotlib.patches import Rectangle

import argparse

# Get stuff for photometry
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
# from photutils import datasets
from photutils.detection import DAOStarFinder
from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.centroids import centroid_sources

# MIRSI data reduction pipeline     Joseph Hora
# Version number of this program
progversion = "v1.30 (2024/02/07)"


# Function for interactive position input in some modules below
def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % 
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))


# Function to determine if a string has real numbers in it
def is_real(item):
    # Method to find if an 'item'is real number
    item = str(item).strip()
    if not(item):
        return False
    elif(item.isdigit()):
        return True
    elif re.match(r"\d+\.*\d*", item) or re.match(r"-\d+\.*\d*", item):
        return True
    else:
        return False


donoise = False
do_ext_calc = False
Jy_ADU = 1.0  # for calibration

argParser = argparse.ArgumentParser(prog='MIRSIphot',
                description='Performs photometry and calibration on MIRSI images',
                epilog='')
argParser.add_argument("-j", "--jyadu", help="Calibration factor (Jy/ADU)",
                       type=float)
argParser.add_argument("-e", "--extfac",
                        help="Extinction factor (override default)")
argParser.add_argument("-s", "--sigma", action='store_true',
                        help="Determine 1 sigma uncertainty from background noise")
argParser.add_argument("-d", "--docal",
                       help="Report flux using these values for Jy, ADU of cal star",
                       nargs=2, type=float)
argParser.add_argument("filename", help="File name of image")

args = argParser.parse_args()

if args.sigma:
    donoise = True
if args.filename:
    filename = args.filename
    if not os.path.exists(filename):
        print("Error: file not found")
        sys.exit(2)
else:
    print("Error - must specify image filename")
    sys.exit(3)
if args.jyadu:
    Jy_ADU = args.jyadu
else:
    Jy_ADU = 1.0
if args.extfac:
    do_ext_calc = False
if args.docal:
    if args.docal[1] > 0:
        Jy_ADU = args.docal[0] / args.docal[1]
        print("\nUsing calibration factor Jy/ADU: ", str("%12.6E\n" % Jy_ADU))
    else:
        print("Error: ADUs and Jy must be greater than zero")
        sys.exit(2)
# Perform photometry on source by clicking on image: Single Image mode

# ADUs are corrected with Airmass in final flux number
# Common                 RA(2000)      Dec(2000)                 Approx. Flux densities (Jy) at various wavelengths (um)
# Name           HD#     (h:m:s)        (d:m:s)       SpType     7.9      8.9      9.8     10.4    11.9    12.9    18.3    Reference
# ------------  ------  ------------  ------------   ----------  -------  -------  ------  ------- ------- ------- ------- ---------
# Beta Cet      4128    00:43:35.372  -17:59:11.77   K0III       85.60    68.69    57.57   49.52   39.44   35.30   16.97   TIMMI2 

# Jy_ADU = 8.406E-5 # Jy/ADU 2/22
# Jy_ADU = 68.69/108235.92 # Jy/ADU for Beta Cet 8.7um
# Jy_ADU = 49.52 /361538.2 # Jy/ADU for Beta Cet 10.57 
# Jy_ADU = 39.44 /63975.963 # Jy/ADU for Beta Cet 11.7 um

# Jy_ADU = 9.613E-5  # Jy/ADU 7/28 alpha Lyr
# Jy_ADU = 1.0855E-4  # Jy/ADU 7/10 alpha Her
# Jy_ADU = 1341/ 12024623  # Jy/ADU 8/04 alpha Her 1.1152116785698812E-4

# Jy_ADU = 169.023/1517073  # for Eta Sgr  (1552211 + 1481935.9)/2   1.1141E-4

# 2022/10/06 observing: NEO 65803
# Jy_ADU = 227/1764906.1  # Beta And 10.57 microns   1.28618E-4
# Jy_ADU = 312/484823.53  # Beta And 8.7    6.435331222E-4
# Jy_ADU = 176/357138.96  # alpha CMa 8.7    4.928053774E-4
# Jy_ADU = 102/202965.65  # alpha CMa 11.7    5.0254809E-4
# Jy_ADU = 119/1202977.8    # alpha CMa 10.57   9.89211e-05

# 20230325 - NEO 2023_   alpha Hya
# Jy_ADU = 115/966078.72     # alpha Hya 10.57 microns  115 Jy (Cohen)

# 20230416 - 
# Jy_ADU = 104.9/995433.95   # beta Gem 10.57 microns (-1.22 Cohen N-band, 104.9 Jy) file 21-40
# Jy_ADU = 104.9/922326.1   # beta Gem 10.57 microns (-1.22 Cohen N-band, 104.9 Jy) file 247-266
# Jy_ADU =  87.9/749526.86  # Mu UMa 10.57 microns (-0.96 10.57 microns 87.9 Jy) file 271-286  Jy/ADU = 0.00011727
# Jy_ADU =  87.9/756194.19  # Mu UMa 10.57 microns (-0.96 10.57 microns 87.9 Jy) file 327-346  Jy/ADU = 0.00011623
# Jy_ADU =  87.9/733218.79  # Mu UMa 10.57 microns (-0.96 10.57 microns 87.9 Jy) file 447-466  Jy/ADU = 0.00011988

# Jy_ADU = (87.9/749526.86  + 87.9/756194.19)/2

# Krisciunas, K., Sinton, W., Tholen, K., Tokunaga, A., Golisch, W., Griep, D., ,
# Journal: Publications of the Astronomical Society of the Pacific, Vol. 99, NO. AUGUST, P. 887, 1987
# summary of extinction values on Mauna Kea 1980-1986
# median airmass corrections:
#     7.8    0.458
#     8.7    0.120
#     9.8    0.151
#     10 N   0.151
#     10.3   0.074
#     11.6   0.081
#     12.5   0.125
#     20     0.419

# am_ext = 0.120    # for 8.7 microns

# am_ext = 0.081    # for 11.7 micron extinction

matplotlib.use('Qt5Agg')

i = 0
x = [1]
y = [1]
wra = [1]
wdec = [1]
imname = ['test']
x.clear()
y.clear()
imname.clear()
wra.clear()
wdec.clear()

image = filename

hdulist = fits.open(image)
hdr = hdulist[0].header
wcs = WCS(hdr)
image_data = hdulist[0].data
airmass_obs = 1.0 #hdr['AMASSAVG']
wl = hdr['LAMBDA']
# Krisciunas, K., Sinton, W., Tholen, K., Tokunaga, A., Golisch, W., Griep, D., ,
# Journal: Publications of the Astronomical Society of the Pacific, Vol. 99, NO. AUGUST, P. 887, 1987
# summary of extinction values on Mauna Kea 1980-1986
# median airmass corrections:
#     7.8    0.458
#     8.7    0.120
#     9.8    0.151
#     10 N   0.151
#     10.3   0.074
#     11.6   0.081
#     12.5   0.125
#     20     0.419

if do_ext_calc:
    if wl == 10.57:
        am_ext = 0.151  # N band extinction
    elif wl == 7.8:
        am_ext = 0.458
    elif wl == 8.7:
        am_ext = 0.120
    elif wl == 9.8:
        am_ext = 0.151
    elif wl == 10.3:
        am_ext = 0.074
    elif wl == 11.6:
        am_ext = 0.081
    elif wl == 12.5:
        am_ext = 0.125
    elif wl > 13:
        am_ext = 0.419

mean, median, std = sigma_clipped_stats(image_data, sigma=3.0)
print((mean, median, std))
image_data = image_data - median

fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(projection=wcs)

cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.title(image)
plt.imshow(image_data, cmap='gray', origin='lower', vmin=-2 * std,
           vmax=5 * std,)
pts = plt.ginput(1)  # number of clicks
print("x:", pts[0][0], pts[0][1])
plt.close()

xref, yref = centroid_sources(image_data, pts[0][0], pts[0][1],
                              box_size=5)
# xref, yref = centroid_sources(image_data, 182, 198, box_size = 9)
positions = (xref[0], yref[0])
print("Positions: ", positions)
apertures = CircularAperture(positions, r=8)
annulus_apertures = CircularAnnulus(positions, r_in=10, r_out=15)
#    apers = [apertures, annulus_apertures]
#    phot_table = aperture_photometry(image_data, apers)

mask = annulus_apertures.to_mask(method='center')

bkg_median = []
annulus_data = mask.multiply(image_data)
annulus_data_1d = annulus_data[mask.data > 0]
_, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
bkg_median.append(median_sigclip)
bkg_median = np.array(bkg_median)
phot = aperture_photometry(image_data, apertures)
phot['annulus_median'] = bkg_median
phot['aper_bkg'] = bkg_median * apertures.area
phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
if do_ext_calc:
    cfactor = Jy_ADU * pow(10, am_ext * airmass_obs / 2.5)
else:
    cfactor = Jy_ADU
phot['flux'] = (phot['aper_sum_bkgsub'] * cfactor)

if do_ext_calc:
    print('Airmass:', str('%5.3f' % airmass_obs))

for col in phot.colnames:
    phot[col].info.format = '%.8g'  # for consistent table output

if Jy_ADU == 1.0:
    print("\nFlux value below is in ADU")
else:
    print("\nFlux value below is in Jy")
print(phot)

w = WCS(hdr)

for ii in range(len(phot)):
    if ((phot['aper_sum_bkgsub'][ii] > 8000)
            and (phot['aper_sum_bkgsub'][ii] < 50000000)):
        y.append(phot['aper_sum_bkgsub'][ii])
        x.append(i)
        wx, wy = w.wcs_pix2world(phot['xcenter'][ii],
                                 phot['ycenter'][ii], 1)
        c = SkyCoord(ra=wx * u.degree, dec=wy * u.degree)
        wra.append(c.to_string('hmsdms', precision=2))
        imname.append(image)
        i = i + 1

# Estimate noise by determining standard deviation of background pixels
if donoise:
    matplotlib.use('Qt5Agg')

    print("\n\nClick twice in different corners to define box")
    print("  click twice at same location to end\n")
    repeat_sample = True

    ref_data = getdata(image)
    ref_hdr = getheader(image)
    wcs_ref = WCS(ref_hdr)
    print("Evaluating image: ", image)
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(projection=wcs)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.title("Select region with two clicks-")
    imstd = np.std(ref_data)
    plt.imshow(ref_data, cmap='gray', vmin=-imstd, vmax=imstd * 4,
               origin='lower')
    while repeat_sample:
        pts = plt.ginput(2, timeout=-1)  # number of clicks
        x1 = pts[0][0]
        y1 = pts[0][1]
        x2 = pts[1][0]
        y2 = pts[1][1]
        # End sampling by clicking on same spot twice
        if abs(pts[0][0] - pts[1][0]) > 3:
            lycorn = int(max(0, min(y1, y2)))
            lxcorn = int(max(0, min(x1, x2)))
            uycorn = int(min(max(y1, y2), ref_data.shape[1]))
            uxcorn = int(min(max(x1, x2), ref_data.shape[0]))
            print("box: ", lxcorn, uxcorn, lycorn, uycorn)
            img_std = np.std(ref_data[lycorn:uycorn, lxcorn:uxcorn])
            plt.gca().add_patch(Rectangle((lxcorn, lycorn),
                                          (uxcorn - lxcorn + 1),
                                          (uycorn - lycorn + 1),
                                          edgecolor='red',
                                          facecolor='none',
                                          lw=1))
            # 1 sigma: per pixel noise times sqrt(number of noise pixels=58.1)
            print("Uncertainty: ",
                  str("%8.5f" % (img_std * Jy_ADU * np.sqrt(58.1))), " Jy\n")
            fig.canvas.draw()
            fig.canvas.flush_events()
        else:
            repeat_sample = False
            print("done.")

    plt.close()
    
