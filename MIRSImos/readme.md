# MIRSImos - make mosaics from sky-subtracted images

The MIRSImos program is run to align the images and create mosaics from the frames. 
The program reads in the subtracted files and takes out some of the array artifacts like column to 
column offsets that look like vertical stripes on the array, and applies an optional flat field 
or "gain" map. The individual corrected frames are written out with a ".cen.fits" at the end of the 
filenames. Then the mosaic is made from these corrected frames. 

The mosaic is made by reprojecting all of the frames to a common WCS, and then averaging the frames with 
sigma clipping to remove bad pixels. The program generates a default output name that includes the object 
name, filter, and range of image files used in the mosaic, and ends with "_mosaic.fits". A second image 
showing the number of overlapping frames is also saved, with the "_coadd.fits" suffix. 
```
> python .\MIRSImos.py  -h
usage: MIRSImos [-h] [-i IRTFCODE] [-d DATECODE] [-o OBJECT] [-p PATH] [-r RANGE RANGE] [-r2 RANGE2 RANGE2] 
                [-m MASK] [-c CENMETHOD] [-n] [-f FILENAME] [-a] [-mr MEDROWS MEDROWS] [-mc MEDCOLS MEDCOLS]
                [-fm FAINTMAX] [-ac] [-s SKYCUT] [-g GAIN]

Makes MIRSI mosaics frames from nod observations

options:
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
  -av --avoidmode       Use the column and row values to define the
                        region to avoid
  -fm FAINTMAX, --faintmax FAINTMAX
                        Faint source maximum value
  -ac, --amasscorr      correct frames for airmass
  -s SKYCUT, --skycut SKYCUT
                        Sky cutoff value (frames rejected if above)
```
The -i, -d, -o, -p, and -r switches control which files are read in for the mosaicing operation. the Gain image can be
used to correct the image for pixel-to-pixel gain variations. This image is multipied with each frame. 

There are three options for aligning the images, set with the -c switch: 

- auto - this works well for bright objects like standard stars which may or may not have had guiding. 
The program uses the first file as a reference image and aligns the other images using 2d cross correlation. 
The object does not need to be a point source. 
- blind - assumes that the TCS offsets in the file headers 
are correct and applies them to the WCS in the image headers. This is the mode to use for the data where one has used MOC
to guide on the source and the header offsets are therefore correct.  
- interactive - each frame is displayed to the screen and the user can 
click on the source in the image to use to align the frames. You click twice on the center of the object
you want to use, and the program calculates a centroid to determine the source position, so you don't have to be exact. 
If a frame is displayed that you don't want to use, because of some problem with the image, you separate 
your clicks by >5 pixels on the image and it will then ignore that frame and not use it in the mosaic. 
The interactive mode is useful if there are multiple sources in the image, or if there are some images that you 
do not want used in the mosaic, or if there are artifacts that are confusing the cross-correlation method and
you need to specify the object in the frame for the program to use for the alignment.

## Column/row offset removal

Another important function of this program is that column and row median offsets are removed before mosaicing. 
By default the program uses almost the entire array, avoiding 5 columns and rows around the edges, to determine the
median level of each column and row. The range can be set with the -mr and -mc switches. For faint point sources, the entire 
array can be used since the few pixels involved with the source will not significantly affect the median in a column or row. 
Another option is to use the -av or avoid mode switch, which will then interpret the range of columns and rows as the region
of the array to avoid using for the median calculation. This is useful for cases where a large extended source is present in the 
middle of the field of view, such as when observing Jupiter. 

## Output file 

By default, the program constructs an output file name based on parameters of the input files. The file is written to
the directory where the input files are located.

The "OBJECT" name in the first FITS file is used for the root of the output file name. If this keyword is not present,
the user is prompted to enter the name. Spaces are replaced by 
underscores to avoid them in the file name, and any parenthesis are removed.
The filter name is appended, along with the range of image numbers used in the mosaic.
The text "_mosaic.fits" is appended to the end. For example, a standard star mosaic name could be "Mu_UMa_10.57_1-20_mosaic.fits".

## Airmass correction
If the data were taken over a short period of time so that the airmass was relatively constant throughout the sequence, such as for
the standard star observations, the user could choose not to apply an airmass correction to the frames and just use the median airmass
for the file set, which is stored in the mosaic header. Then the airmass correction could be performed in the photometry step. However,
if the airmass does change significantly during the observation, e.g. if the object is setting at low altitude or one just has a long sequence 
of frames, then an airmass correction can be applied to the individual frames before mosaicing. The program uses extinction values based on
Krisciunas, K., Sinton, W., Tholen et al. 1987 (PASP, 99, 887), which are as follows:
```
            median airmass
     filter corrections
----------------------------
     7.8    0.458
     8.7    0.120
     9.8    0.151
     10 N   0.151
     10.3   0.074
     11.6   0.081
     12.5   0.125
     20     0.419
```

## Mosaic construction
The program uses the `find_optimal_celestial_wcs` and `reproject_exact` procedures from the `reproject` package to reproject each image
to a common WCS, and then combines the images using averaging with sigma clipping (sigma=2). A "coadd.fits" image is also saved which records the
number of overlapping frames at each point of the output mosaic. 
