# MIRSIdiff - make sky-subtracted MIRSI images

This program is the first step in reducing MIRSI data. It performs the A-B nod frame subtraction which removes most of the sky background.
The MIRSI files are of the format:

```
mrs.YYYYSNNN.YYMMDD.object.nnnnn.a.fits
```

where "mrs." is the prefix given to all MIRSI data files, YYYYSNNN is the program ID giving the year, semester (A or B) and the 
program number NNN, "object" is the text for the file name entered by the observer 
in the MIRSI observing program, "nnnnn" is the file number, and the end of the file is either "a.fits" or "b.fits" for the A or B beam. 
For example, "mrs.2024A029.240320.setup.00001.a.fits" is a file obtained in program 2024A029 on 2024/03/20.

MIRSIdiff subtracts the images and writes both files, with the B file being the inverse of the A file.
The output file names are the same as the input files but with a ".sub.fits" at the end. Below are the command line options:

```
> python MIRSIdiff.py  -h
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
The -i, -d, and -o switches are for setting those values in the file name for selection of the files to operate on. The range of file 
numbers can be defined by the -r switch. If -l is selected, a log file of the observations is created. 

```
MIRSI Observing log                        Date: 2023-03-19
                                        Program: 2023A072
                                      Observers: Trilling,Lopez-Oquendo,Hora

TIME(UT)   FRNUM     OBJECT     FILESTR  FILTER  COADDS   ITIME  BEAM   RAOFF   DEOFF   AIRMASS
---------------------------------------------------------------------------------------------------
05:43:18     1     beta Gem         std   10.57     200   0.015     a   0.00    0.00     1.023
05:43:28     2     beta Gem         std   10.57     200   0.015     b  -6.00  -15.00     1.023
05:43:38     3     beta Gem         std   10.57     200   0.015     a -10.40    3.30     1.022
05:43:48     4     beta Gem         std   10.57     200   0.015     b -16.40  -11.70     1.022
05:43:59     5     beta Gem         std   10.57     200   0.015     a  -3.10   -2.80     1.022
05:44:08     6     beta Gem         std   10.57     200   0.015     b  -9.10  -17.80     1.022
05:44:19     7     beta Gem         std   10.57     200   0.015     a -13.40   -2.30     1.022
05:44:29     8     beta Gem         std   10.57     200   0.015     b -19.40  -17.30     1.022
05:44:39     9     beta Gem         std   10.57     200   0.015     a  -0.80    3.50     1.022
05:44:49    10     beta Gem         std   10.57     200   0.015     b  -6.80  -11.50     1.022
```
