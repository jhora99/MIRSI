
# # MIRSImos - MIRSI image reduction software

# Import all necessary packages and functions

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import re
import sys

# from astropy.utils.data import download_file
from astropy.io import fits
from astropy.io.fits import getval
from astropy.io.fits import getheader
from astropy.wcs import WCS


# photometry /image alignment functions
from photutils.centroids import centroid_sources
from image_registration import chi2_shift

# For mosaic construction
# from reproject import reproject_interp
from reproject import reproject_exact
# from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs
from astropy.stats import sigma_clip
import warnings
warnings.filterwarnings('ignore')

import argparse

# MIRSI data reduction pipeline     Joseph Hora
# Version number of this program
progversion = "v1.41 (2024/02/14)"

# Version history:
#  v1.1 (2022/01/19) - Initial version
#        - offsets are read from COMMENT line in data headers (put in by observing macro)
#        - beamswitches hardwired in, not put in headers by MIRSI observing program
#        - using plate scale information from original MIRSI calibration
#        - now normalizing to ADU/s in the initial A-B sky subtraction section
#
#  v1.15 (2022/02/03) - first NEO run
#        - sign of offsets was changed
#  v1.16 (2022/02/22) - second NEO run
#        - sign of offsets was changed
#  v1.20 (2022/04/11) - Jupiter runs
#        - added option to avoid a range of pixels in the column-wise and row-wise medians
#

# Function for interactive position input in some modules below
def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))


# Function to determine if a string has real numbers in it

def is_real(item):
#""" Method to find if an 'item'is real number"""
    item = str(item).strip()
    if not(item):
        return False
    elif(item.isdigit()):
        return True
    elif re.match(r"\d+\.*\d*", item) or re.match(r"-\d+\.*\d*", item):
        return True
    else:
        return False

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
def get_extfactor(wl):
    if wl>6 and wl<8:
        cval = 0.458
    elif wl>=8 and wl<9:
        cval = 0.120
    elif wl==10.57:
        cval = 0.151
    elif wl>9 and wl<10:
        cval = 0.151
    elif wl>=10 and wl<11:
        cval = 0.074
    elif wl>=11 and wl<12:
        cval = 0.081
    elif wl>=12 and wl<13:
        cval = 0.125
    elif wl>16 and wl<21:
        cval = 0.419
    else:
        print('Warning - extinction not defined, setting to 0.08')
        cval = 0.08
    return cval

def amass_factor(wl, amass):
    return pow(10, get_extfactor(wl) * amass / 2.5)
       
    

# This function returns the filter name and wavelength of the filter wheel setup
# for the broadband filters, this is mostly from wheel 1, but there are a couple filters in wheel 2
# This does not check for non-standard setups, for example 2.2 in w1 and blank in w2
# would return the 2.2 filter and wavelength, even though there would be no transmission
def get_filter(wheel1, wheel2):
    wh1filters = ["Q0", "Q1", "2.2", "20.7", "Q2", "12.282", "10.57", "12.5", "11.7", "9.8", "7.7", "4.9", "8.7"]
    wv1filters = [16.8, 17.5,   2.2,  20.7,   18.4, 12.282,   10.57,   12.5,   11.7,   9.8,   7.9,   4.9,   8.7]
    wh2filters = ["Q3", "K", "24.8", "Blank", "20uGrism", "10uGrism",]
    wv2filters = [19,   2.2, 24.8,     0,          20,       10]
    if wheel2 == "Blank":
        wavelength = 0.0
        filtername = "Blank"
    elif wheel1 in wh1filters:
        filtername = wheel1
        wavelength = wv1filters[wh1filters.index(filtername)]
    elif wheel2 in wh2filters:
        filtername = wheel2
        wavelength = wv2filters[wh2filters.index(filtername)]
    return filtername, wavelength
    

# This function contains the possible processes that can be done to a sky-subtracted 
# MIRSI frame to remove medians, flip the data, apply gain, zero out data outside of
# a defined range. Inputs are the image and header plus the flags and parameters
# that specify the operations to be done. The header is updated by adding HISTORY
# lines as a record of the operations performed.

def process_mirsiframe(image_data, hdr, subflip, submed, doGain, ColMedSub, RowMedSub, 
                       NegMax, PosMax, YCmin, YCmax, XRmin, XRmax, AvoidMode):
    if subflip:
        image_data = -image_data
        hdr.add_history('IMAGE DATA WAS MULTIPLIED BY (-1)')
    if submed:
        medval = np.median(image_data[YCmin:YCmax,XRmin:XRmax])
        image_data = np.subtract(image_data,medval)
        hdr.add_history('IMAGE MEDIAN SUBTRACTED FROM IMAGE: ' + str("%8.3f" % medval))
    
    if ColMedSub:
        if AvoidMode:   # in this mode, do not use pixels in the region defined by YCmin to YCmax
            calc_array = image_data[1:YCmin,:]
            calc_array2 = image_data[YCmax+1:240,:]
            calc_array = np.concatenate((calc_array, calc_array2))
            med_columns = np.median(calc_array, axis=0)
            image_data = image_data - med_columns
            hdr.add_history('COLUMN MEDIANS SUBTRACTED, region avoided:' + 
                            str("%3i" % YCmin) + ', ' + str("%3i" % YCmax))
        else:           # in this mode, use pixels in the range YCmin to YCmax
            calc_array = image_data[YCmin:YCmax,:]
            med_columns = np.median(calc_array, axis=0)
            image_data = image_data - med_columns
            hdr.add_history('COLUMN MEDIANS SUBTRACTED: min and max used:' + 
                            str("%3i" % YCmin) + ' to ' + str("%3i" % YCmax))
    if RowMedSub:
        if AvoidMode:      # in this mode, do not use pixels in the region defined by XRmin to XRmax
            image_t = image_data.T
            calc_array = image_t[1:XRmin,:]
            calc_array2 = image_t[XRmax:320,:]
            calc_array = np.concatenate((calc_array, calc_array2))
            med_rows = np.median(calc_array, axis=0)
            calc_array = image_t - med_rows
            image_data = calc_array.T    
            hdr.add_history('ROW MEDIANS SUBTRACTED, region avoided:' + 
                            str("%3i" % XRmin) + ' to ' + str("%3i" % XRmax))
        else:             # in this mode, use pixels in the range XRmin to XRmax
            image_t = image_data.T
            calc_array = image_t[XRmin:XRmax,:]
            med_rows = np.median(calc_array, axis=0)
            calc_array = image_t - med_rows
            image_data = calc_array.T    
            hdr.add_history('ROW MEDIANS SUBTRACTED: min and max used:' + 
                            str("%3i" % XRmin) + ', ' + str("%3i" % XRmax))
    if doGain:
        image_data = image_data * gain_image
        hdr.add_history('GAIN IMAGE USED: ' + gain_name)

    image_data = np.where(image_data>NegMax, image_data, np.nan)
    image_data = np.where(image_data<PosMax, image_data, np.nan)
    hdr.add_history('IMAGE VALUES SET TO NAN OUTSIDE RANGE ' + 
                    str("%6.0f" % NegMax) + " to " + str("%6.0f" % PosMax))
    return image_data, hdr


# Define the location and filenames to search for processing

propcode = "2023A072"
datecode = "230416"
objcode = "*"              # use * for any name after the datecode in the filename
impath = "./"
cenmethod = 'blind'
irange2 = False 
NoMosaic = False
amass_corr = False
faintmode = False
# Set the options below to control which frames are processed and how

# Select the filenames and (optionally) the filter setting to process
filter_select = ''     # leave blank if no selection on filter needed
registration_text = "" # this gets filled in later with registration method

# Select by file number
# You can restrict the range of the file numbers to use
# set inum_max to zero if no check required
inum_min = 0       # minimum file number in sequence
inum_max = 5000    # maximum file number in sequence

# If you want to add a second range of file numbers to process this dataset, set irange2=True
irange2 = False      # use a second range of images?  Set True if wanted
inum_min2 = 0        # minimum file number in sequence
inum_max2 = 0        # maximum file number in sequence
# inum_max = inum_min + 11   
# The following flag will make the limits specified above to be the region to AVOID
# This is useful when looking at large extended objects like Jupiter that nearly fill the frame

# If true, calculate medians for rows and columns everywhere EXCEPT in values above
AvoidRegion = False     

# Set True if source is in both A and B frames. B frame will be
# multiplied by -1 and used in the mosaic
AandB = True
# multiply differenced frames by -1: this was necessary for early versions of
# instrument control program
subflip = False         

# Masking pixels by their value in the background-subtracted frames
#smallest valid data value in subtracted frame, mask values below this
NegMax =  -10000
# largest allowed value in subtracted frame; mask values above this
PosMax = 1591000


# Controls how the images are processed 
RowMedSub = True       # subtract row-wise median
ColMedSub = True       # subtract column-wise median
FrameMedSub = True     # subtract frame-wise median


YCmin = 20             # Min and max Y values to use in column median calculation
YCmax = 200            # (max<240)
XRmin = 90             # Min and max columns 
XRmax = 250            # over which to calculate row medians (max<320)

doMask = False
doGain = False
do_skycut = False
useroutname = ""

# Read in command line arguments 
argParser = argparse.ArgumentParser(prog='MIRSImos',
                    description='Makes MIRSI mosaics frames from nod observations',
                    epilog='')
argParser.add_argument("-i", "--irtfcode", help="IRTF Proposal code")
argParser.add_argument("-d", "--datecode", help="Observing date code (YYMMDD)")
argParser.add_argument("-o", "--object", help="Object string in filename (default is *)")
argParser.add_argument("-p", "--path", help="path to data files")
argParser.add_argument("-r", "--range", 
                       help="Range of file numbers to process", nargs=2, type=int)
argParser.add_argument("-r2", "--range2", 
                       help="Optional second range to use", nargs=2, type=int)
argParser.add_argument("-g", "--gain", help="Gain image file name")
argParser.add_argument("-m", "--mask", 
                       help="Value in gain image below which to mask pixels", type=float)
argParser.add_argument("-c", "--cenmethod", 
                       help="Centering method (auto, blind, or interactive)")
argParser.add_argument("-n", "--nomosaic", 
                       help="Do not make mosaic of frames", action='store_true')
argParser.add_argument("-f", "--filename",
                       help="name of output mosaic (overrides default name)")
argParser.add_argument("-a", "--Aframe", help="A frames only", action='store_true')
argParser.add_argument("-mr", "--medrows", help="Range for Median calculations for rows", nargs=2, type=int)
argParser.add_argument("-mc", "--medcols", help="Range for Median calculations for columns", nargs=2, type=int)
argParser.add_argument("-fm", "--faintmax", help="Faint source maximum value", type=float)
argParser.add_argument("-ac", "--amasscorr", help="correct frames for airmass", action='store_true')
argParser.add_argument("-s", "--skycut", help="Sky cutoff value (frames rejected if above)", type=float)

args = argParser.parse_args()

if args.amasscorr:
    amass_corr = True
if args.faintmax:
    faintmode = True
    faintmax = args.faintmax
if args.irtfcode:
    propcode = args.irtfcode
if args.datecode:
    datecode = args.datecode
if args.object:
    objcode = args.object
if args.path:
    impath = args.path
if args.range:
    inum_min = min(args.range)
    inum_max = max(args.range)
if args.range2:
    irange2 = True
    inum_min2 = min(args.range2)
    inum_max2 = max(args.range2)
if args.gain:
    doGain = True          # use gain map
    gain_name = args.gain
if args.mask:
    doMask = True          # mask bad pixels when making mosaics
    maskval = args.mask
if args.cenmethod:
    cenmethod = args.cenmethod
if args.nomosaic:
    NoMosaic = True
if args.filename:
    useroutname = args.filename
if args.Aframe:
    AandB = False
if args.medcols:
    YCmin = max(min(args.medcols),0)             # Min and max Y values to use in column median calculation
    YCmax = min(max(args.medcols),240)           # (max<240)
if args.medrows:
    XRmin = max(min(args.medrows),0)             # Min and max columns 
    XRmax = min(max(args.medrows),319)           # over which to calculate row medians (max<320)
if args.skycut:
    skycutoff = args.skycut
    do_skycut = True

    
if doGain:
    hdulist=fits.open(gain_name)
    gain_image = hdulist[0].data
    # flip gain image up/down and left/right if necessary
    gain_image = np.fliplr(np.flipud(gain_image))
    hdulist.close()
if doMask:
    # construct mask map from low response pixels
    mask_image = np.where(gain_image>maskval,1,0)    
    mask_name = "GAIN <" + str("%8.5f" % maskval)


# The following section constructs a list of frames to process, given the specifications above. 
# The files in the imdir are searched and a list of files is produced in mirsilist

# If the source is in both nod beams, include the ".b." frames from the list in the reduction
if AandB:
    matchlist = glob.glob(impath + "/mrs." + propcode + "." + datecode + 
                          "." + objcode + '*.?.sub.fits')
else:
    matchlist = glob.glob(impath + "/mrs." + propcode + "." + datecode + 
                          "." + objcode + '*.a.sub.fits')

# If the filter_select flag is set, only add those frames that match the desired filter
if len(filter_select)>0:
    mirsilist = []
    for image in matchlist:
        
        if getval(image, "FILTER", 0) == filter_select:
            mirsilist.append(image)
else:
    mirsilist = matchlist.copy()

testlist = mirsilist.copy()

if len(mirsilist)==0:
    print("Warning - no files found!")
    sys.exit(2)
    
# set "final" path character to either a forward or backward slash, depending on OS
testname = mirsilist[0]

if testname.rfind("/") > testname.rfind("\\"):
    fchar = "/"
else:
    fchar = "\\"
    
# Remove all frames from the list that are outside the requested range (no removals if inum_max=0)
if inum_max>0:
    if irange2:
        for image in testlist:
            inum = int(image[len(image)-16:len(image)-11])
            if not ((inum <= inum_max and inum >= inum_min) or 
                    (inum <= inum_max2 and inum >= inum_min2)):
                mirsilist.remove(image)
    else:
        for image in testlist:
            inum = int(image[len(image)-16:len(image)-11])
            if not (inum <= inum_max and inum >= inum_min):
                mirsilist.remove(image)

if len(mirsilist)==0:
    print("Warning - no files found!")
    sys.exit(2)

print("Object name: ", getval(mirsilist[0],'OBJECT'))
# show the list of files selected
# print(mirsilist)


# # Below are three possible options for aligning the frames. Use one of these before 
# making the mosaics

# Option 1 =  cross-correlate the images: first frame is the reference, all other 
#             frames aligned to the first one

# Option 2 =  use blind offsetting based on the telescope offsets (good for when 
#             guiding with MOC)

# Option 3 =  interactive, click on the position of the source to align, and a local
#              centroid is calculated to determine offset

# Initialize value for summing the airmasses (in order to average them later)
am_sum = 0.0
# sum for calculating average observation time
mjd_sum = 0.0

if (cenmethod == 'auto'):
    #  OPTION 1:
    
    # Automatic offsets and image processing by registering the frames based on the reference
    # image using an image cross-correlation technique
    #
    # This is useful for standard star frames with no guiding. The star is typically the brightest
    # source in the frame and this routine will work well to determine the relative offsets
    #
    # The images are cross-correlated with the reference image to determine the offset.
    # This is then used to fix the reference position in the file header.
    #
    # The RA, Dec in the first frame in the sequence is used as the reference position.
    #
    # All images are assumed valid and used in the subsequent mosaic generation.
    print("Using cross-correlation alignment")
    goodlist = []
    ref_image = mirsilist[0]
    print("Reference image: ",ref_image)
    registration_text = "CROSS CORRELATION TO REFERENCE IMAGE USED TO ALIGN FRAMES"
    
    refhdulist = fits.open(ref_image, readonly = True)
    ref_data = refhdulist[0].data
    
    ref_data, refhdulist[0].header = process_mirsiframe(ref_data, refhdulist[0].header, subflip, FrameMedSub, 
                                                           doGain, ColMedSub, RowMedSub, NegMax, PosMax, 
                                                           YCmin, YCmax, XRmin, XRmax, AvoidRegion)
    
    ref_x = refhdulist[0].header['CRPIX1'] 
    ref_y = refhdulist[0].header['CRPIX2'] 
    raref = refhdulist[0].header['CRVAL1'] 
    decref = refhdulist[0].header['CRVAL2'] 
    ref_hdr = refhdulist[0].header
    
    for objimname in mirsilist:
        hdulist = fits.open(objimname, readonly = True)
        hdr = hdulist[0].header
        image_data = hdulist[0].data
        image_data, hdulist[0].header = process_mirsiframe(image_data, hdulist[0].header, subflip, FrameMedSub, 
                                                           doGain, ColMedSub, RowMedSub, NegMax, PosMax, 
                                                           YCmin, YCmax, XRmin, XRmax, AvoidRegion)
    
        x1, y1, exoff, eyoff = chi2_shift(ref_data, image_data, upsample_factor='auto')
        hdulist[0].data = image_data
        hdulist[0].header['CRPIX1'] = ref_x + x1
        hdulist[0].header['CRPIX2'] = ref_y + y1
        hdulist[0].header['CRVAL1'] = raref
        hdulist[0].header['CRVAL2'] = decref
        am_sum += hdulist[0].header['TCS_AM']
        mjd_sum += hdulist[0].header['MJD_OBS']
        outimage = objimname[:objimname.find('sub.fits')]+'sub.cen.fits'
        print(outimage[outimage.find('mrs'):], 'offsets:', str("%5.2f" % x1),str("%5.2f" % y1))
        hdulist.writeto(outimage, overwrite=True)
        if do_skycut:
            if hdulist[0].header['I_MEDIAN']<skycutoff and hdulist[0].header['S_MEDIAN']<skycutoff:
                goodlist.append(outimage)
        else:   
            goodlist.append(outimage)

        hdulist.close()

elif cenmethod == 'blind':
    # OPTION 2:
    
    # Register images with "blind offsets" and image processing using header coordinates - 
    #  This routine is useful when using MOC guiding, for moving or stationary object
    #  The offsets are determined from the commanded position as saved in the header
    #
    #  For early data taking:
    #  The offsets are read from the header comments, the beamswitch must be hardwired
    #  in for moving objects. This is no longer needed as of 2022 observing runs. 
    #
    # All images are assumed valid and used in the subsequent mosaic generation.
    print("Using blind offset image alignment")
    goodlist = []
    registration_text = "HEADER OFFSETS AND BEAMSWITCH VALUES USED TO ALIGN FRAMES"
    for objimname in mirsilist:
        hdulist = fits.open(objimname, readonly = True)
        hdr = hdulist[0].header
        if objimname == mirsilist[0]:
            refhdr = hdulist[0].header
        if objimname == mirsilist[1]:
            skyhdr = hdulist[0].header
        image_data = hdulist[0].data
        image_data, hdulist[0].header = process_mirsiframe(image_data, hdulist[0].header, subflip, FrameMedSub, 
                                                           doGain, ColMedSub, RowMedSub, NegMax, PosMax, 
                                                           YCmin, YCmax, XRmin, XRmax, AvoidRegion)
        x1 = hdulist[0].header['CRPIX1']
        y1 = hdulist[0].header['CRPIX2']
        commenttext = hdulist[0].header['COMMENT']
        offlist=[]
        xoff = 0
        yoff = 0
        
    
        # New method - read total offsets from header, convert to pixels    
        if hdulist[0].header.get('OS_TRA', default=-9999) != -9999:
            # Convert arcsec offsets to units of pixels
            xoff = hdulist[0].header['OS_TRA']/0.269
            yoff = hdulist[0].header['OS_TDEC']/0.264
        else:
            print('Warning: Total offsets not present in MIRSI headers!')
        hdulist[0].header['CRVAL1'] = refhdr['CRVAL1']
        hdulist[0].header['CRVAL2'] = refhdr['CRVAL2']
        total_yoff = refhdr['CRPIX2'] - yoff
        total_xoff = refhdr['CRPIX1'] + xoff
        hdulist[0].header['CRPIX1'] = xoff
        hdulist[0].header['CRPIX2'] = -yoff
        am_sum += hdulist[0].header['TCS_AM']
        mjd_sum += hdulist[0].header['MJD_OBS']
        outimage = objimname[:objimname.find('sub.fits')]+'sub.cen.fits'
        hdulist[0].data = image_data
        print(outimage[outimage.find('mrs'):], 'offsets:', 
              str("%5.2f" % xoff),str("%5.2f" % yoff), str("%5.2f" % total_xoff),
              str("%5.2f" % total_yoff))
        hdulist.writeto(outimage, overwrite=True)
        if do_skycut:
            print('image med: ',str('%9.1f' % hdulist[0].header['I_MEDIAN']), 
                  'sky med: ', str('%9.1f' % hdulist[0].header['S_MEDIAN']), end="")
            if hdulist[0].header['I_MEDIAN']<skycutoff and hdulist[0].header['S_MEDIAN']<skycutoff:
                goodlist.append(outimage)
                print(" ")
            else:
                print(" - rejected")
        else:   
            goodlist.append(outimage)
        hdulist.close()

elif cenmethod == 'interactive':
    # OPTION 3:
    
    # Interactive offsets and image processing
    #
    # The processed images are displayed one at a time, and the user clicks on the 
    # position of the source. Two clicks in nearly the same position mark the object 
    # and accept it for use in the mosaic. Two clicks at positions >5 pixels
    # different will mark the frame as not to be used.
    #
    # Then a centroid is calculated in a box around the marked position to get a precise 
    # offset. This is then used to fix the reference position in the file header.
    #
    # The RA, Dec in the first frame in the sequence is used as the reference position.
    print("Interactive image alignment:\n")
    matplotlib.use('Qt5Agg')
    goodlist=[]
    registration_text = "FRAMES ALIGNED INTERACTIVELY BY USER"
    print('Click on source twice to use for alignment')
    print('If clicks are >5 pixels apart, frame is rejected')
    for ref_image in mirsilist:
        refhdulist = fits.open(ref_image, readonly = True)
        ref_data = refhdulist[0].data
        ref_data, refhdulist[0].header = process_mirsiframe(ref_data, refhdulist[0].header, subflip, FrameMedSub, 
                                                           doGain, ColMedSub, RowMedSub, NegMax, PosMax, 
                                                           YCmin, YCmax, XRmin, XRmax, AvoidRegion)
        if ref_image == mirsilist[0]:
            ref_hdr = refhdulist[0].header
            raref = refhdulist[0].header['CRVAL1'] 
            decref = refhdulist[0].header['CRVAL2'] 
            ref_x = refhdulist[0].header['CRPIX1'] 
            ref_y = refhdulist[0].header['CRPIX2'] 
        wcs = WCS(refhdulist[0].header)
        print("Evaluating image: " + ref_image[ref_image.rfind(fchar)+1:])
        fig = plt.figure(figsize=(12,12))
        ax=plt.subplot(projection=wcs)
        
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.title(ref_image)
        imstd = np.nanstd(ref_data)
        plt.imshow(ref_data, cmap='gray', vmin=-imstd, vmax=imstd*4, origin='lower') 
        pts = plt.ginput(2) #number of clicks
        print("position selected:     " + str("%6.3f" % pts[0][0]) + ", " 
              + str("%6.3f" % pts[0][1]))
        
        # If the two selected points are close, use the image
        if abs(pts[0][0]-pts[1][0])<5:
            xref, yref = centroid_sources(ref_data, pts[0][0],  pts[0][1], box_size = 15)
            print("Offsets from centroid: " + str("%6.3f" % xref[0]) + "," 
                  + str("%6.3f" % yref[0]))
            refhdulist[0].header['CRPIX1'] = xref[0]
            refhdulist[0].header['CRPIX2'] = yref[0]
            refhdulist[0].header['CRVAL1'] = raref 
            refhdulist[0].header['CRVAL2'] = decref
            refhdulist[0].data = ref_data
            outimage = ref_image[:ref_image.find('sub.fits')]+'sub.cen.fits'
            if do_skycut:
                if hdulist[0].header['I_MEDIAN']<skycutoff and hdulist[0].header['S_MEDIAN']<skycutoff:
                    goodlist.append(outimage)
            else:   
                goodlist.append(outimage)
            am_sum += refhdulist[0].header['TCS_AM']
            mjd_sum += hdulist[0].header['MJD_OBS']
            refhdulist.writeto(outimage, overwrite=True)
            refhdulist.close()
        else:
            print(ref_image[ref_image.rfind(fchar)+1:] + " rejected, not used")
        plt.close()
else:
    print('Invalid method')
    sys.exit(3)


#****************************************************************************************
#
# Make Mosaic from frames
#
#****************************************************************************************

# From the list of images in "goodlist", combine the images based on the WCS in the files

if NoMosaic:
    print("no mosaic mode")
elif len(goodlist) > 0:
    print("now making mosaic....")
    
    wcs_out, shape_out = find_optimal_celestial_wcs(goodlist)
    print("Mosaic size:", shape_out[1],"x",shape_out[0])
    image_concat = []
    fprint_concat = []
    refimnum = max(0,int(len(goodlist)/2)-1)
    ref_hdr = getheader(goodlist[refimnum])
    ref_name = goodlist[refimnum][(goodlist[refimnum].rfind(fchar)+1):]
    
    # Get the filter, object name, and coadds used to construct the output file name
    filtername, wavelength = get_filter(ref_hdr['GFLT'], ref_hdr['CVF'])
    objname = ref_hdr['OBJECT']
    coadds = ref_hdr['CO_ADDS']
    
    for image in goodlist:
        hdu_list = fits.open(image)
        sourceframe = hdu_list[0].data
        header = hdu_list[0].header
        amass_val = header['TCS_AM']
        if faintmode or amass_corr:
            sourceframe = sourceframe - np.nanmedian(sourceframe)
        if faintmode:
            sourceframe = np.where(np.abs(sourceframe)<faintmax, sourceframe, np.NAN)
        if amass_corr:
            sourceframe = sourceframe * amass_factor(wavelength, amass_val)
        if doMask:
            sourceframe = np.where(mask_image==0,np.nan,sourceframe)
        hdu_list[0].data = sourceframe
        array, footprint = reproject_exact(hdu_list, wcs_out, shape_out=shape_out)
        image_concat.append(array)
        np.seterr(invalid='ignore')
        farray = np.divide(array, array)        # make array with 1=good, nan for bad pixels
        fprint_concat.append(farray)
    
    stacked_array = np.ma.stack(image_concat)
    stacked_fprint = np.ma.stack(fprint_concat)
    coadd_image = np.nansum(stacked_fprint, axis=0)
    filtered_data = sigma_clip(stacked_array, axis=0, sigma=2, masked=False, copy=False)#, maxiters=None,
                                #cenfunc=mean)#, masked=False, copy=False)
    med_image = np.nanmean(filtered_data, axis=0)
    med_image = np.where(coadd_image<1, np.nan, med_image)
    
    header = wcs_out.to_header()
    
    # zero out NaNs
    med_image = np.where(np.isnan(med_image), 0, med_image)
    
    hdum = fits.PrimaryHDU(med_image, header=header)
    hdul = fits.HDUList([hdum])
    
    # Copy header information from first file in list to the output mosaic header so that it 
    # will retain useful information like filter, object name, comments, etc.
    
    # Must first check to see if header keywords are present already in new header (for example
    # the WCS info), otherwise those cards will be overwritten
    
    for fitskeywd in ref_hdr:
        found = False
        for wcsfitskeywd in hdul[0].header:
            if wcsfitskeywd == fitskeywd:
                found = True
        if not found:                      # If the keyword is not already in header, add it
            if fitskeywd == "COMMENT":                # The comment and history keywords have 
                commenttext = ref_hdr[fitskeywd]      # multi-line values, so must add them
                for line in commenttext:              # line-by-line  
                    hdul[0].header.add_comment(line)
            elif fitskeywd == "HISTORY":
                commenttext = ref_hdr[fitskeywd]
                for line in commenttext:
                    hdul[0].header.add_history(line)
            else:
                hdul[0].header.append((fitskeywd, ref_hdr[fitskeywd], 
                                       ref_hdr.comments[fitskeywd]), end=True)
    
    # Add additional header information specifying how mosaic was made
    if doMask:
        hdul[0].header.add_history("IMAGES MASKED WITH " + mask_name)
    hdul[0].header.add_history("IMAGE NUMBER RANGE: " 
                               + str("%i"% inum_min) + " - " + str("%i"% inum_max))
    if irange2:
        hdul[0].header.add_history("SECOND IMAGE NUMBER RANGE: " 
                                   + str("%i"% inum_min2) + " - " + str("%i"% inum_max2))
    
    hdul[0].header.add_history("HEADER INFORMATION ABOVE FROM MIDDLE FILE OF SET:")
    hdul[0].header.add_history("REF NAME:" + ref_name)
    if len(registration_text)>0:
        hdul[0].header.add_history(registration_text)
    hdul[0].header.append(('NMOSFRAM', len(goodlist),'Number of frames used to make mosaic'))
    hdul[0].header.add_history("IMAGES WERE COMBINED WITH SIGMA-CLIPPED MEAN (SIGMA=2)")
    if amass_corr:
        hdul[0].header.add_history("AIRMASS CORRECTION APPLIED TO INDIVIDUAL FRAMES")
        hdul[0].header['AMASSEXT'] = (get_extfactor(wavelength), 'EXTINCTION VALUE USED IN AMASSCORR')
    if faintmode:
        hdul[0].header.add_history("MEDIANS SUBTRACTED FROM INDIVIDUAL FRAMES")
        hdul[0].header.add_history(("DATA VALUES CLIPPED, ABS(DATA)<"+str('%10.5E' % faintmax)))
    hdul[0].header['AMASSAVG'] = (am_sum / len(goodlist), 'AVERAGE AIRMASS IN MOSAIC')
    hdul[0].header['MJD_OBS'] = (mjd_sum / len(goodlist), 'AVERAGE MJD OF OBSERVATION')
    outfname = objname.replace(" ","_")

    # Make file name from object, filter, and frame numbers with Ncoadds
    outfile = impath + '/' + outfname # + '_' + str(coadds) 
    outfile = outfile + '_' + filtername + "_" + str(inum_min) + '-' + str(inum_max)
    if irange2:
        outfile = outfile + "_" + str(inum_min2) + '-' + str(inum_max2)
    # remove parentheses and spaces from filename if present
    outfile = outfile.replace('(', '_')
    outfile = outfile.replace(')', '_')
    outfile = outfile.replace(' ', '_')
    outfile = outfile.replace('__','_')
    hdul = fits.HDUList([hdum])
    if len(useroutname) > 0:
        outfile = useroutname
    # write out the mosaic file
    print("writing mosaic file: ", (outfile + '_mosaic.fits'))
    hdul.writeto(outfile + '_mosaic.fits', overwrite=True)
    
    # # write out coadd file
    hdum = fits.PrimaryHDU(coadd_image.data, header=header)
    hdul = fits.HDUList([hdum])
    hdul.writeto(outfile + '_coadd.fits', overwrite=True)
else:
    print("Error: No images in mosaic list")
