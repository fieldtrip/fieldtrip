Souheil J. Inati
Dartmouth College
May 2000
souheil.inati@dartmouth.edu

Listing of GE2SPM directory:

GE_convertADW.m       GE_dumpSeriesHeader.m  GE_readHeaderSuite.m
GE_convertVolume.m    GE_readHeader.m        GE_readVolume.m
GE_createSPMHeader.m  GE_readHeaderExam.m    GE_reorientHeader.m
GE_dumpExamHeader.m   GE_readHeaderImage.m   GE_reorientImage.m
GE_dumpHeader.m       GE_readHeaderPixel.m   GE_writeSPMHeader.m
GE_dumpImageHeader.m  GE_readHeaderSeries.m  README

The main scripts are GE_convertVolume, GE_convertADW, and GE_dumpHeader.

The conversion scripts reshape and flip the data around to match what
SPM expects, ie axial in radiological convention.  This has been
tested with Ax,Sag,Cor with slices acquired SI and IS, LR and RL, and
PA, PA.  Also tested for Oblique axial.  Don't count on double
obliques or anything really fancy

GE_convertVolume
  Converts a series of 5.X GE slices into Analyze format for SPM.
Should be easy to convert this to handle LX images.

GE_convertADW
  Converts a series of GE slices acquired on the Advanced Development
Workstation into Analyze format for SPM.  These are 5.X format stored
on the Sun machine GE uses to bypass the 512 images/series limitation
of its database.

GE_dumpHeader
 Dumps much of the header from a GE lx2 or 5.X file to a file or stdout.
