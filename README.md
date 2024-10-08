# MIRSI reduction programs
The Mid-Infrared Spectrometer and Imager (MIRSI) is an instrument in use at the [NASA Infrared Telescope Facility](https://irtfweb.ifa.hawaii.edu/) (IRTF). It can obtain images in the 2-20 micron range, and optical images at 0.5 - 1 micron with an integrated CCD camera that views the same field through a dichroic mirror.

See the [MIRSI paper](https://ui.adsabs.harvard.edu/abs/2024arXiv240902752H) for a description of MIRSI and some example images obtained with the camera.

See the [MIRSI reduction cookbook](https://github.com/jhora99/MIRSI/blob/main/MIRSIcookbook.md) for an 
example of how to reduce MIRSI data.

See the readme files in the program directories for more information on each module.

## [MIRSIdiff](./MIRSIdiff)
Create a log of observations, and/or perform the A-B nod subtraction operation to remove sky background

## [MIRSImos](./MIRSImos)
Combine MIRSI sky-subtracted frames into a mosaic

## [MIRSIphot](./MIRSIphot)
Perform aperture photometry on MIRSI images

## [data](./data)
The data directory contains the transmission or reflection curves for various filters and optical elements in MIRSI.
