# # MIRSIphot - MIRSI image reduction and photometry software

# MIRSI data reduction pipeline     Joseph Hora

# Import all necessary packages and functions

import numpy as np
import glob

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io.fits import getval

import argparse

# Version number of this program
progversion = "v1.4 (2024/05/10)"

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
#  v1.30 (2023/04/17) convert from jupyter notebook to python
#
#  v1.4 (2024/05/10) - add log file construction


# This function returns the filter name and wavelength of the filter wheel setup
# for the broadband filters, this is mostly from wheel 1, but there are a couple
# filters in wheel 2
# This does not check for non-standard setups, for example 2.2 in w1 and blank
# in w2 would return the 2.2 filter and wavelength, even though there would be
# no transmission

def get_filter(wheel1, wheel2):
    wh1filters = ["Q0", "Q1", "2.2", "20.7", "Q2", "12.282", "10.57",
                  "12.5", "11.7", "9.8", "7.7", "4.9", "8.7"]
    wv1filters = [16.8, 17.5,   2.2,  20.7,   18.4, 12.282,   10.57,
                  12.5,   11.7,   9.8,   7.9,   4.9,   8.7]
    wh2filters = ["Q3", "K", "24.8", "Blank", "20uGrism", "10uGrism"]
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
    
    
propcode = "*"
datecode = "*"
objcode = "*"              # use * for any name after the datecode in the filename
impath = '.'
inum_max = 0
verbose = False
logfile = False
writelog=False

# Define the location and filenames to search for processing

argParser = argparse.ArgumentParser(prog='MIRSIdiff',
                description='Makes MIRSI difference frames from nod observations',
                epilog='')
argParser.add_argument("-i", "--irtfcode", help="IRTF Proposal code")
argParser.add_argument("-d", "--datecode", help="Observing date code (YYMMDD)")
argParser.add_argument("-o", "--object",
                       help="Object string in filename (default is *)")
argParser.add_argument("-p", "--path", help="path to data files")
argParser.add_argument("-r", "--range",
                       help="Range of file numbers to process",
                       nargs=2, type=int)
argParser.add_argument("-v", "--verbose",
                       help="Print out list of files while processing",
                       action='store_true')
argParser.add_argument("-l", "--logfile", 
                       help="Write MIRSI obseration log file (no frame processing)")

args = argParser.parse_args()

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
if args.verbose:
    verbose = True
if args.logfile:
    writelog = True
    logname = args.logfile
    
imdir = impath + "/"
image_root = imdir + "mrs." + propcode + "." + datecode + "." + objcode + "."
searchpath = image_root+'0????.?.fits'
print('Looking for MIRSI files in ', searchpath)

# Construct a list of files to perform A-B sky subtraction

# The first code block at the end of this section assumes that frames are
# taken in AB mode, and the sky subtraction will be frame A - frame B. The
# B - A frame is also saved and if the source is in both beams, the B-A
# frame can be used in the mosaic by setting flags in cells below
#
# If the data are taken as consecutive frames in beam A, then use the other
# code block at the end of this section to set up the lists of frames
# to subtract

imlist = []
skylist = []
itimes = []
objects = []
coadds = []
filters = []
inums = []
beams = []
filelist = []
filenamelist = []
times = []
raoff = []
deoff = []
airmasses = []
medflux = []
allfilelist = glob.glob(searchpath)


def last_12chars(x):
    return(x[-12:])


# sort filenames by file number to get into correct order for sky subtraction
allfilelist = sorted(allfilelist, key=last_12chars)
fnstr = '.' + datecode + '.'

# Look at the files in list and get object names, itime, coadds, filter, beam,
# frame number
for image in allfilelist:
    hdu_list = fits.open(image, readonly=True)
    filtername, wavelength = get_filter(hdu_list[0].header['GFLT'],
                                        hdu_list[0].header['CVF'])
    beammode = hdu_list[0].header['BEAMPAT']
    normfactor = hdu_list[0].header['CO_ADDS']*hdu_list[0].header['ITIME']
    medflux.append(np.median(hdu_list[0].data)/normfactor)
    if beammode > 0:               # make list of files for sky subtraction
        filters.append(filtername)             # with only files in AB mode
        objects.append(hdu_list[0].header['OBJECT'])
        itimes.append(hdu_list[0].header['ITIME'])
        coadds.append(hdu_list[0].header['CO_ADDS'])
        inums.append(int(image[len(image) - 12:len(image) - 7]))
        beams.append(image[len(image) - 6:len(image) - 5])
        times.append(hdu_list[0].header['TIME_OBS'][:8])
        airmasses.append(hdu_list[0].header['TCS_AM'])
        if hdu_list[0].header.get('OS_TRA')==None:
            raoff.append(0)
            deoff.append(0)   
        else:    
            raoff.append(hdu_list[0].header['OS_TRA'])
            deoff.append(hdu_list[0].header['OS_TDEC'])
        filelist.append(image)
        tmpstr = image[:-13]
        filenamelist.append(tmpstr[(tmpstr.rfind('.')+1):])
    hdu_list.close()

# Make two lists, one with A frames (object) and the other with B frames
# (sky or object)
print('Total number of frames:     ', len(filelist))

if verbose:
    print('Frnum  object  filenm  filter coadds itime  beam filename')
if writelog:
    # Open log file and write header
    logfile = open(impath + '/' + logname, "w")
    logfile.write("MIRSI Observing log                        Date: " +
                  getval(filelist[0], 'DATE_PC'))
    logfile.write("\n                                        Program: " +
                  getval(filelist[0],'PROG_ID'))
    logfile.write("\n                                      Observers: " +
                  getval(filelist[0],'OBSERVER'))
    logfile.write("\n\nTIME(UT)   FRNUM     OBJECT     FILESTR  FILTER  COADDS   ITIME  BEAM   RAOFF   DEOFF   AIRMASS    MEDIAN\n")
    logfile.write("-------------------------------------------------------------------------------------------------------------\n")
# Frames assumed to be taken in AB mode (alternating beams A and B images)
for idx in range(0, len(filenamelist)):
    if verbose:
        print(inums[idx], objects[idx], filenamelist[idx], filters[idx],
              coadds[idx], itimes[idx], beams[idx], filelist[idx])
    if writelog:
        logline = times[idx] + str("%6i " % inums[idx]) + str("%13s" % objects[idx])
        logline = logline + str("%12s" % filenamelist[idx])
        logline = logline +  str("%8s" % filters[idx]) + str("%8i" % coadds[idx])
        logline = logline + str("%8.3f" % itimes[idx]) + str("%6s" % beams[idx])
        logline = logline + str("%7.2f " % raoff[idx]) + str("%7.2f" % deoff[idx])
        logline = logline + str("%10.3f " % airmasses[idx]) + str("%10.1f" % medflux[idx]) + '\n'
        logfile.write(logline)
    if beams[idx] == 'a':
        imlist.append(filelist[idx])
    else:
        skylist.append(filelist[idx])

        
# If writing the observing log, close file and end the program. Otherwise,
# proceed to do the frame differences that were specified

if writelog:
    logfile.close()
else:    
    # If a range was specified, remove file names outside of that range
    if inum_max > 0:
        testlist = imlist.copy()
        for image in testlist:
            inum = int(image[len(image) - 12:len(image) - 7])
            if not (inum <= inum_max and inum >= inum_min):
                imlist.remove(image)
        testlist = skylist.copy()
        for image in testlist:
            inum = int(image[len(image) - 12:len(image) - 7])
            if not (inum <= inum_max and inum >= inum_min):
                skylist.remove(image)

    # Find the type of "slash" character at the end of the path: PCs 
    # might have "\", otherwise it is "/"
    if len(imlist) > 0:
        if imlist[0].rfind("/") > imlist[0].rfind("\\"):
            fchar = "/"
        else:
            fchar = "\\"


    # SKY SUBTRACTION

    # This does the A-B differencing of the frames. It assumes that there are
    # pairs for all frames, and the closest sky frame is subtracted. If doing
    # 1A-2B-3A-4B-5A-6B etc. pattern, it does
    #   (1A-2B), (3A-4B), (5A-6B)
    # if the data are in 1A-2B-3B-4A-5A-6B etc. pattern, it would do
    #   (1A-2B), (4A-3B), (5A-6B)

    # The data are also normalized by a factor of (itime * coadds) to make the
    # images ADU/sec
    # Some header information is fixed/added
    i = 0
    for image, sky in zip(imlist, skylist):
        hdusky = fits.open(sky, readonly=True)
        hdulist = fits.open(image, readonly=True)
        hdr = hdulist[0].header
        hdrsky = hdusky[0].header
        image_data = hdulist[0].data
        image_median = np.nanmedian(image_data)
        skydata = hdusky[0].data
        sky_median = np.nanmedian(skydata)
        imnum = int(image[len(image)-12:len(image)-7])
        skynum = int(sky[len(sky)-12:len(sky)-7])
        if abs(imnum - skynum) > 1:
            print("WARNING: non-consecutive frame subtraction!")
        # subtract background from object frame
        diffarray = np.subtract(image_data, skydata)
        outfile = image[:image.find("fits")]+'sub.fits'
        outsky = sky[:sky.find("fits")]+'sub.fits'
        # flip data up/down to orient it properly
        diffarray = np.flipud(diffarray)
        normfactor = hdr['CO_ADDS']*hdr['ITIME']
        diffarray = diffarray / normfactor               # normalize to ADU/sec
        image_median = image_median / normfactor
        sky_median = sky_median / normfactor
        # remove first 20 columns for bad readout board
        diffarray = diffarray[:,20:]       
        mfilter, wavelength = get_filter(hdr['GFLT'], hdr['CVF'])
        if verbose:
            print(outfile[outfile.rfind(fchar)+1:] + " - "
                  + outsky[outsky.rfind(fchar)+1:])
        # put WCS info into proper header keywords
        hdu = fits.PrimaryHDU(diffarray, header=hdr)
        headcoord = hdu.header['TCS_RA'] + hdu.header['TCS_DEC']
        c = SkyCoord(headcoord, unit=(u.hourangle, u.deg))
        hdu.header.append(('FILTER', mfilter, 'Broadband filter name'))
        hdu.header.append(('LAMBDA', wavelength, 'Wavelength of filter'))
        hdu.header.append(('CTYPE1', 'RA---TAN', 'WCS of image'))
        hdu.header.append(('CTYPE2', 'DEC--TAN'))
        hdu.header.append(('CRVAL1', c.ra.degree))
        hdu.header.append(('CRVAL2', c.dec.degree))
        hdu.header.append(('CROTA2', 0.0))
        hdu.header.append(('CDELT1', -7.472E-05))
        hdu.header.append(('CDELT2', 7.333E-05))
        hdu.header.append(('CRPIX1', 160))
        hdu.header.append(('CRPIX2', 120))
        hdu.header.append(('EQUINOX', 2000.0))
        hdu.header.append(('I_MEDIAN', image_median,'IMAGE MEDIAN'))
        hdu.header.append(('S_MEDIAN', sky_median,'SKY MEDIAN'))
        hdu.header.add_history("MIRSI PIPELINE VERSION: " + progversion)
        hdu.header.add_history("THE FOLLOWING IMAGE WAS SUBTRACTED TO REMOVE SKY BACKGROUND")
        hdu.header.add_history("FILE: " + outsky[outsky.rfind(fchar)+1:])
        hdul = fits.HDUList([hdu])
        hdul.writeto(outfile, overwrite=True)  # write the sky-subtracted image
        hdul.close()
        # Write out the B beam frame as well, in case the source is there too
        # but use the inverse of the subtraction so B beam source will be positive
        hdu = fits.PrimaryHDU(-diffarray, header=hdrsky)
        headcoord = hdu.header['TCS_RA'] + hdu.header['TCS_DEC']
        c = SkyCoord(headcoord, unit=(u.hourangle, u.deg))
    #     hdu.header['INSTRUME'] = 'MIRSI'
        hdu.header.append(('FILTER', mfilter, 'Broadband filter name'))
        hdu.header.append(('LAMBDA', wavelength, 'Wavelength of filter'))
        hdu.header.append(('CTYPE1', 'RA---TAN', 'WCS of image'))
        hdu.header.append(('CTYPE2', 'DEC--TAN'))
        hdu.header.append(('CRVAL1', c.ra.degree))
        hdu.header.append(('CRVAL2', c.dec.degree))
        hdu.header.append(('CROTA2', 0.0))
        hdu.header.append(('CDELT1', -7.472E-05))
        hdu.header.append(('CDELT2', 7.333E-05))
        hdu.header.append(('CRPIX1', 160))
        hdu.header.append(('CRPIX2', 120))
        hdu.header.append(('EQUINOX', 2000.0))
        hdu.header.append(('I_MEDIAN', image_median,'IMAGE MEDIAN'))
        hdu.header.append(('S_MEDIAN', sky_median,'SKY MEDIAN'))
        hdu.header.add_history("MIRSI PIPELINE VERSION: " + progversion)
        hdu.header.add_history("THE FOLLOWING IMAGE WAS SUBTRACTED TO REMOVE SKY BACKGROUND")
        hdu.header.add_history("FILE: " + outfile[outfile.rfind(fchar)+1:])
        hdul = fits.HDUList([hdu])
        hdul.writeto(outsky, overwrite=True)    # write the sky-subtracted image
        hdul.close()

print('Processing complete.')
