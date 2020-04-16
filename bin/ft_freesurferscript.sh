#!/bin/sh

# This script sets up the freesurfer environment and executes 
# freesurfer's automatic processing pipeline, resulting in a
# set of meshes that model the cerebral cortex. It takes two
# input arguments. Use as:
#
#  ft_freesurferscript.sh <DIRNAME> <SUBJECTNAME> 
# 
#  1. <DIRNAME> is the directory that contains a <SUBJECTNAME>.mgz
#      anatomical image, that will be used as input volume
#  2. The <SUBJECTNAME>. The script produces a subfolder in <DIRNAME>.
#      called <SUBJECTNAME>, which contains the results of the processing
#
# The script can be ran without user intervention, but it may require
# manual intervention at specific stages, due to the automatic steps 
# yielding suboptimal results. Most often this is caused by suboptimal
# white matter segmentation. Users may want to inspect the outcome of the
# volumetric part of the processing pipeline (-autorecon1 and part of 
# -autorecon2) before proceeding to the remainder (rest of -autorecon2 and
# -autorecon3). Please refer to the freesurfer documentation for more info
# about this. Note that the script works for freesurfer 6.0, it may need 
# some tweaks in order for it to work for other versions.

# Copyright (C), 2019, Jan-Mathijs Schoffelen

# update the directory of the freesurfer code if needed
export FREESURFER_HOME=/opt/freesurfer/6.0
export SUBJECTS_DIR=$1

source $FREESURFER_HOME/SetUpFreeSurfer.sh

mksubjdirs $SUBJECTS_DIR/$2

# This step copies the T1w anatomical image
cd $SUBJECTS_DIR
cp -f $2.mgz $SUBJECTS_DIR/$2/mri/

# This steps converts the image as per the freesurfer-required
# input format
cd $SUBJECTS_DIR/$2/mri
mri_convert -c -oc 0 0 0 $2.mgz orig.mgz
cp orig.mgz orig/001.mgz

# Run
recon-all -autorecon1 -subjid $2
recon-all -autorecon2 -subjid $2
recon-all -autorecon3 -subjid $2
