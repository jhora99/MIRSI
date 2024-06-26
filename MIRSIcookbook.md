
# Cookbook for MIRSI data reduction                         
### Joseph Hora
### 2023/04/20

This document shows some examples of how to reduce MIRSI data for the NEO observing program. 

The first program 20230418_MIRSIdiff.py performs the A-B subtraction, to give you basically 
what you see on the screen when observing. It writes out both files, the B file is the inverse 
of the A file, to be used to make the mosaic or do photometry on to get a positive flux. The 
files are the same name as the input file but with a "sub.fits" at the end. 

The second program 20230418_MIRSImos.py reads in the subtracted files and takes out some of 
the array artifacts like column to column offsets that look like vertical stripes on the 
array, and applies an optional flat field or "gain" map. I have attached an example gain map 
that I use that we made from some calibration observations taken during engineering time. 
The program also aligns the images using one of three different modes, and makes a mosaic 
of the frames. The individual corrected frames are written out with a ".cen.fits" at the 
end of the filenames. 

The third program 20230419_MIRSIphot.py performs photometry on point sources in the images. 
You first run this on standard stars to determine the Jy/ADU factor. Then you can run the 
program on the NEO frames to get their flux. The program also does the airmass correction 
based on the average airmass read from the mosaic header (calculated by the MIRSImos program 
in the previous step). The program optionally measures the noise in the image and outputs an estimate of the 1-sigma point source sensitivity level.

If you run the programs with the -h option they will print out the information on what command line switches and options you use to specify the files to process and various modes of the programs. 


## Step 1. Perform the image differencing and save the results as "*.sub.fits" files 


This step uses the MIRSIdiff python program. Here is an example of the processing for the 230207 data:

```
> python .\20230418_MIRSIdiff.py  -h
usage: MIRSIdiff [-h] [-i IRTFCODE] [-d DATECODE] [-o OBJECT] [-p PATH] [-r RANGE RANGE] [-v] [-l LOGFILE]

Makes MIRSI difference frames from nod observations

optional arguments:
  -h, --help            show this help message and exit
  -i IRTFCODE, --irtfcode IRTFCODE
                        IRTF Proposal code
  -d DATECODE, --datecode DATECODE
                        Observing date code (YYMMDD)
  -o OBJECT, --object OBJECT
                        Object string in filename (default is *)
  -p PATH, --path PATH  path to data files
  -r RANGE RANGE, --range RANGE RANGE
                        Range of file numbers to process
  -v, --verbose         Print out list of files while processing
  -l LOGFILE, --logfile LOGFILE
                        Write MIRSI obseration log file
```

For this dataset, there was apparently an issue around file number 58-59 where two beam A frames were recorded back-to-back. 
In order for the A-B frames to be subtracted correctly, the data were processed in the following two steps. If all of the 
data were in A,B,A,B, etc. pattern, then the entire night could have been done with one command without specifying a range.

```
> python .\20230418_MIRSIdiff.py -i 2023A072 -d 230207 -p e:\MIRSI\NEOobs\230207 -v -l logfile.txt -r 1 58
> python .\20230418_MIRSIdiff.py -i 2023A072 -d 230207 -p e:\MIRSI\NEOobs\230207 -v -r 60 166
```

Using the "-v" option, the file names and results are printed to the screen. The "-l logfile.txt" option made a log file 
for the entire night's data (that is not affected by the range option). The log file shows useful information including 
file numbers, object names, airmass, times, filters, exposure times, etc.


## Step 2: Removing artifacts and aligning images, making mosaics

Then the MIRSImos program is run to align the images and create mosaics from the frames. The program reads in the subtracted files and takes out some of the array artifacts like column to column offsets that look like vertical stripes on the array, and applies an optional flat field or "gain" map. The gain map that I use ("gain.1.2.fits") was made from some calibration observations taken during engineering time.  The individual corrected frames are written out with a ".cen.fits" at the end of the filenames. Then the mosaic is made from these corrected frames. 

There are three options for aligning the images. One is the "auto" method, which works well for bright objects like standard stars which may or may not have had guiding. The program uses the first file as a reference image and aligns the other images using 2d cross correlation. The object does not need to be a point source. The second option for aligning images is the "blind" mode, where it assumes that the TCS offsets in the file headers are correct and applies them to the WCS in the image headers. This is the mode to use for the NEO data where we are  guiding with MOC so we know the header offsets are correct.  

The third alignment mode is "interactive", where each frame is displayed to the screen and the user can click on the source in the image to use to align the frames. You click twice on the center of the object you want to use, and the program calculates a centroid to determine the source position, so you don't have to be exact. If a frame is displayed that you don't want to use, because of some problem with the image, you separate your clicks by >5 pixels on the image and it will then ignore that frame and not use it in the mosaic. The interactive mode is useful if there are multiple sources in the image, or if there are some images that you do not want used in the mosaic, or if there are artifacts that are confusing the cross-correlation method and you need to specify the object in the frame for the program to use for the alignment.

The mosaic is made by reprojecting all of the frames to a common WCS, and then averaging the frames with sigma clipping to remove bad pixels. The program generates a default output name that includes the object name, filter, and range of image files used in the mosaic, and ends with "_mosaic.fits". A second image showing the number of overlapping frames is also saved, with the "_coadd.fits" suffix. See the [MIRSImos](MIRSImos/readme.md) readme file for more information.
```
> python .\20230418_MIRSImos.py  -h
usage: MIRSImos [-h] [-i IRTFCODE] [-d DATECODE] [-o OBJECT] [-p PATH] [-r RANGE RANGE] [-r2 RANGE2 RANGE2] [-g GAIN]
                [-m MASK] [-c CENMETHOD] [-n] [-f FILENAME] [-a] [-mr MEDROWS MEDROWS] [-mc MEDCOLS MEDCOLS]
                [-av AVOIDMODE] [-fm FAINTMAX] [-ac] [-s SKYCUT]

Makes MIRSI mosaics frames from nod observations

optional arguments:
  -h, --help            show this help message and exit
  -i IRTFCODE, --irtfcode IRTFCODE
                        IRTF Proposal code
  -d DATECODE, --datecode DATECODE
                        Observing date code (YYMMDD)
   -o OBJECT, --object OBJECT
                        Object string in filename (default is *)
  -p PATH, --path PATH  path to data files
  -r RANGE RANGE, --range RANGE RANGE
                        Range of file numbers to process
  -r2 RANGE2 RANGE2, --range2 RANGE2 RANGE2
                        Optional second range to use
  -g GAIN, --gain GAIN  Gain image file name
  -m MASK, --mask MASK  Value in gain image below which to mask pixels
  -c CENMETHOD, --cenmethod CENMETHOD
                        Centering method (auto, blind, or interactive)
  -n, --nomosaic        Do not make mosaic of frames
  -f FILENAME, --filename FILENAME
                        name of output mosaic (overrides default name)
  -a, --Aframe          A frames only
  -mr MEDROWS MEDROWS, --medrows MEDROWS MEDROWS
                        Range for Median calculations for rows
  -mc MEDCOLS MEDCOLS, --medcols MEDCOLS MEDCOLS
                        Range for Median calculations for columns
  -av AVOIDMODE, --avoidmode AVOIDMODE
                        Avoid the region defined by the columns and row settings
  -fm FAINTMAX, --faintmax FAINTMAX
                        Faint source maximum value
  -ac, --amasscorr      correct frames for airmass
  -s SKYCUT, --skycut SKYCUT
                        Sky cutoff value (frames rejected if above)
```

Here are the commands used to process two standard stars and one of the sets of NEO data from the 230207 night's data:

```
> python .\20230418_MIRSImos.py -i 2023A072 -d 230207 -p e:\MIRSI\NEOobs\230207 -r 1 20 -c auto
Object name:  Mu_UMa
Using cross-correlation alignment
Reference image:  e:\MIRSI\NEOobs\230207\mrs.2023A072.230207.setup.00001.a.sub.fits
mrs.2023A072.230207.setup.00001.a.sub.cen.fits offsets:  0.00 -0.00
mrs.2023A072.230207.setup.00002.b.sub.cen.fits offsets: -21.90 57.18
mrs.2023A072.230207.setup.00003.a.sub.cen.fits offsets: -39.41 -12.24
.
.
.
mrs.2023A072.230207.setup.00018.b.sub.cen.fits offsets: -13.60 46.91
mrs.2023A072.230207.setup.00019.a.sub.cen.fits offsets: -47.63  0.44
mrs.2023A072.230207.setup.00020.b.sub.cen.fits offsets: -69.50 58.46
now making mosaic....
Mosaic size: 415 x 324
writing mosaic file:  e:\MIRSI\NEOobs\230207/Mu_UMa_10.57_1-20_mosaic.fits


> python .\20230418_MIRSImos.py -i 2023A072 -d 230207 -p e:\MIRSI\NEOobs\230207 -r 130 149 -c auto
Reference image:  e:\MIRSI\NEOobs\230207\mrs.2023A072.230207.setup.00130.a.sub.fits
mrs.2023A072.230207.setup.00130.a.sub.cen.fits offsets:  0.00 -0.00
mrs.2023A072.230207.setup.00131.b.sub.cen.fits offsets: -22.43 57.48
mrs.2023A072.230207.setup.00132.a.sub.cen.fits offsets: -40.09 -12.42
.
.
.
mrs.2023A072.230207.setup.00147.b.sub.cen.fits offsets: -13.43 46.22
mrs.2023A072.230207.setup.00148.a.sub.cen.fits offsets: -47.01  0.19
mrs.2023A072.230207.setup.00149.b.sub.cen.fits offsets: -69.62 58.41
now making mosaic....
Mosaic size: 414 x 324
writing mosaic file:  e:\MIRSI\NEOobs\230207/Mu_UMa_10.57_130-149_mosaic.fits

> python .\20230418_MIRSImos.py  -i 2023A072 -d 230207 -p e:\MIRSI\NEOobs\230207 -r 21 58 -c blind
Object name:  98943
Using blind offset image alignment
mrs.2023A072.230207.neo.00021.a.sub.cen.fits offsets:  0.00  0.00 160.00 120.00
mrs.2023A072.230207.neo.00022.b.sub.cen.fits offsets: -22.30 -56.82 -22.30 56.82
mrs.2023A072.230207.neo.00023.a.sub.cen.fits offsets: -38.66 12.50 -38.66 -12.50
mrs.2023A072.230207.neo.00024.b.sub.cen.fits offsets: -60.97 -44.32 -60.97 44.32
.
.
.
mrs.2023A072.230207.neo.00056.b.sub.cen.fits offsets: -56.51 -54.55 -56.51 54.55
mrs.2023A072.230207.neo.00057.a.sub.cen.fits offsets: 11.15 11.74 11.15 -11.74
mrs.2023A072.230207.neo.00058.b.sub.cen.fits offsets: -11.15 -45.08 -11.15 45.08
now making mosaic....
Mosaic size: 414 x 323
writing mosaic file:  e:\MIRSI\NEOobs\230207/98943_10.57_21-58_mosaic.fits
```
### Faintmax value
If the source is known or expected to have a maximum flux less than a certain ADU level, the "faintmax" value can be 
set with the -fm switch to flag as bad any pixel in the image that is above a certain cutoff level. This is useful for eliminating
noise spikes in the image that can be sometimes much larger than the source, but not rejected by the sigma clipping and so make 
their way into the final image as spurious peaks. This filtering is done for each frame separately, so a pixel marked as 
bad in one image is not necessarily bad in other images. 

### SKYCUT value
If the sky background level is changing throughout the observations, one might not want to use frames where the background
level is above a certain cutoff. This is useful for nights where an occasional thin cloud will pass in front of the telescope
and is visible by an increase in the sky background level. The raw frame median level is checked, and if it is above this 
cutoff value set with the -s switch, the frame is excluded from the mosaic calculation.

### Airmass correction
The airmass correction can be done on individual files in the mosaicing process, or in the photometry steps below. If the object 
has been observed over a wide range of airmasses, it is best to correct it in the mosaicing step (using the -ac option). If
the observation was short, e.g. for a standard star, and the airmass was almost constant during the observation, the correction
could be done in the photometry step. If done in the mosaicing process, a note is made in the header and the AIRMASS value in 
the header is set to zero. If no correction is done while mosaicing, the program determines the average airmass during the 
observation and sets the AIRMASS in the header to the average value. The photometry program will then use this value when performing the correction in later steps.


## Step 3: Photometry of the MIRSI mosaics

The MIRSIphot.py program performs photometry on the mosaics. Any photometry program could be used for this purpose, but this program has
some convenient features for mid-IR photometry.

First run the program on the standard star mosaic(s) to measure the flux in ADUs. Then use that number, along with the known flux of the star in Jy at that wavelength to determine the Jy/ADU factor. Then that factor can be applied to the photometry for the NEO. The airmass correction is done as part of the photometry and is included in the numbers reported.

```
> python .\20230419_MIRSIphot.py  -h
usage: MIRSIphot [-h] [-j JYADU] [-e EXTFAC] [-s] [-d DOCAL DOCAL] filename

Performs photometry and calibration on MIRSI images

positional arguments:
  filename              File name of image

optional arguments:
  -h, --help            show this help message and exit
  -j JYADU, --jyadu JYADU
                        Calibration factor (Jy/ADU)
  -e EXTFAC, --extfac EXTFAC
                        Extinction factor (override default)
  -s, --sigma           Determine 1 sigma uncertainty from background noise
  -d DOCAL DOCAL, --docal DOCAL DOCAL
                        Report flux using these values for Jy, ADU of cal star
```

First, do the photometry on a standard star. When the program is started, a window pops up showing the mosaic image. Click on the object you want to photometer, and the flux in ADUs is returned.

```
> python .\20230419_MIRSIphot.py  e:/MIRSI/NEOobs/230207/Mu_UMa_10.57_1-20_mosaic.fits
(-0.7750801363442089, 0.0, 68.44108638185821)
single click: button=1, x=445, y=547, xdata=170.854839, ydata=189.345161
x: 170.85483870967738 189.34516129032255
Positions:  (170.92896239823514, 190.38206681207961)
Airmass: 1.134

Flux value below is in ADU
 id  xcenter   ycenter  aperture_sum ...  aper_bkg aper_sum_bkgsub    flux
       pix       pix                 ...
--- --------- --------- ------------ ... --------- --------------- ---------
  1 170.92896 190.38207    730781.59 ... 5932.0532       724849.54 848675.11
```

If you want, you can re-run the program on the standard star using the known flux in Jy and the ADU value measured, and you can see that the flux of the star is returned in Jy.

```
> python .\20230419_MIRSIphot.py -d 87.9 848675.11 e:/MIRSI/NEOobs/230207/Mu_UMa_10.57_1-20_mosaic.fits

Using calibration factor Jy/ADU:  1.035732E-04

(-0.7750801363442089, 0.0, 68.44108638185821)
single click: button=1, x=443, y=546, xdata=169.783871, ydata=188.809677
x: 169.7838709677419 188.80967741935478
Positions:  (170.66738655059666, 190.3918857029861)
Airmass: 1.134

Flux value below is in Jy
 id  xcenter   ycenter  aperture_sum ...  aper_bkg aper_sum_bkgsub    flux
       pix       pix                 ...
--- --------- --------- ------------ ... --------- --------------- ---------
  1 170.66739 190.39189    730672.01 ... 5883.9738       724788.03 87.892542
```
Below is a sample image of a standard star:
![MIRSI image of calibration star](f1.png)

Using the "-d" option allows you to enter the flux of the object in Jy and the ADU value returned by the first step. Alternatively, you could instead use the "-j" option which allows you to directly enter the Jy/ADU factor to use.

Now run the program with the same values, but on the NEO mosaic, and you will get the flux of the NEO in Jy. In this example, I have clicked on some random location because I don't see the NEO. It returns a low flux value. Also in this example, I have included the "-s" option, which will calculate the 1-sigma point source sensitivity value in the image. You do this by clicking twice to specify boxes in the image that have no point source, and it determines the sensitivity based on the noise in that section of the image. You can do this for many parts of the image to see that it is consistent in regions where you would expect the NEO to be. Clicking twice at the same location ends the process.

```
> python .\20230419_MIRSIphot.py -d 87.9 848675.11 -s e:/MIRSI/NEOobs/230207/98943_10.57_60-129_mosaic.fits

Using calibration factor Jy/ADU:  1.035732E-04

(-0.31854646004589027, 0.0, 34.61636598280449)
single click: button=1, x=463, y=507, xdata=180.057419, ydata=167.410323
x: 180.05741935483871 167.4103225806452
Positions:  (182.32071397740282, 165.67311555811713)
Airmass: 2.041

Flux value below is in Jy
 id  xcenter   ycenter  aperture_sum ...  aper_bkg aper_sum_bkgsub     flux
       pix       pix                 ...
--- --------- --------- ------------ ... --------- --------------- -----------
  1 182.32071 165.67312    642.19432 ... 365.50627       276.68805 0.038063981


Click twice in different corners to define box
  click twice at same location to end

Evaluating image:  e:/MIRSI/NEOobs/230207/98943_10.57_60-129_mosaic.fits
single click: button=1, x=521, y=557, xdata=211.040645, ydata=194.120000
single click: button=1, x=555, y=513, xdata=229.203226, ydata=170.615484
box:  211 229 170 194
Uncertainty:   0.01469  Jy

single click: button=1, x=409, y=514, xdata=151.210968, ydata=171.149677
single click: button=1, x=475, y=457, xdata=186.467742, ydata=140.700645
box:  151 186 140 171
Uncertainty:   0.01481  Jy

single click: button=1, x=434, y=569, xdata=164.565806, ydata=200.530323
single click: button=1, x=485, y=532, xdata=191.809677, ydata=180.765161
box:  164 191 180 200
Uncertainty:   0.01523  Jy

single click: button=1, x=485, y=532, xdata=191.809677, ydata=180.765161
double click: button=1, x=485, y=532, xdata=191.809677, ydata=180.765161
done.
```
The image below shows an example of background regions selected by clicking in opposite corners of a box
that you want to define:
![MIRSI image of calibration star](f2.png)
Running the photometry on another set of Mu UMa calibration data, it gives a different flux for the star. This indicates that the calibration has changed between the two measurements, probably due to sky conditions. 


```
> python .\20230419_MIRSIphot.py -d 87.9 848675.11 e:/MIRSI/NEOobs/230207/Mu_UMa_10.57_130-149_mosaic.fits

Using calibration factor Jy/ADU:  1.035732E-04

(0.39606424859037137, 0.0, 60.75842771067333)
single click: button=1, x=462, y=559, xdata=179.523226, ydata=195.188387
x: 179.5232258064516 195.18838709677425
Positions:  (180.51288937086727, 194.99553433169328)
Airmass: 1.342

Flux value below is in Jy
 id  xcenter   ycenter  aperture_sum ...  aper_bkg aper_sum_bkgsub    flux
       pix       pix                 ...
--- --------- --------- ------------ ... --------- --------------- ---------
  1 180.51289 194.99553    583456.78 ... 6598.6927       576858.09 72.006748

```
