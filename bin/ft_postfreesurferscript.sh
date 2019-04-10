#!/bin/bash

# This script is adjusted from HCP's FS2CaretConvertRegisterNonlinear
# script, and creates sets of cortical meshes at various resolutions,
# which are surface-registered to a template, and which can be directly
# used for group-level averaging. The results will be stored in the 
# directory 'workbench', in the subject-specific freesurfer folder.
# 
# Use as:
#
#  ft_postfreesurferscript.sh <DIRNAME> <SUBJECTNAME> <TEMPLATEDIRNAME>
#
# 1. <DIRNAME> points to the directory which contains the subject specific
#     directory that contains the freesurfer results.
# 2. <SUBJECTNAME> is the name of the subject specific directory.
# 3. <TEMPLATEDIRNAME> is the directory that contains the required template
#     meshes.
#
# The template directory should be set up in advance by the user, and should
# consist of the Pipelines/global/templates/standard_mesh_atlases directory
# that can be obtained from https://github.com/Washington-University/HCPpipelines.
# In addition, one should add the templates for the low-resolution spherical meshes
# from fieldtrip/template/sourcemodel, i.e. all L/R.*.gii files.
# In addition, this script requires HCP-workbench to be installed.

# Copyright (C), 2019, Jan-Mathijs Schoffelen

set -e
echo -e "\n START: FS2CaretConvertRegisterNonlinear"

Subject="$2"
FreeSurferFolder="$1"/"$Subject"
SurfaceAtlasDIR="$3"

LowResMeshes=4@8@32
LowResMeshes=`echo ${LowResMeshes} | sed 's/@/ /g'`
GrayordinatesResolutions=`echo ${GrayordinatesResolutions} | sed 's/@/ /g'`

WBFolder="$FreeSurferFolder"/workbench
HighResMesh=164

#Make some folders for this and later scripts
if [ ! -e "$WBFolder" ] ; then
  mkdir -p "$WBFolder"
fi
if [ ! -e "$WBFolder"/fsaverage ] ; then
  mkdir -p "$WBFolder"/fsaverage
fi

#Find c_ras offset between FreeSurfer surface and volume and generate matrix to transform surfaces
MatrixX=`mri_info "$FreeSurferFolder"/mri/brain.finalsurfs.mgz | grep "c_r" | cut -d "=" -f 5 | sed s/" "/""/g`
MatrixY=`mri_info "$FreeSurferFolder"/mri/brain.finalsurfs.mgz | grep "c_a" | cut -d "=" -f 5 | sed s/" "/""/g`
MatrixZ=`mri_info "$FreeSurferFolder"/mri/brain.finalsurfs.mgz | grep "c_s" | cut -d "=" -f 5 | sed s/" "/""/g`
echo "1 0 0 ""$MatrixX" > "$FreeSurferFolder"/mri/c_ras.mat
echo "0 1 0 ""$MatrixY" >> "$FreeSurferFolder"/mri/c_ras.mat
echo "0 0 1 ""$MatrixZ" >> "$FreeSurferFolder"/mri/c_ras.mat
echo "0 0 0 1" >> "$FreeSurferFolder"/mri/c_ras.mat

#Loop through left and right hemispheres
for Hemisphere in L R ; do
  #Set a bunch of different ways of saying left and right
  if [ $Hemisphere = "L" ] ; then 
    hemisphere="l"
    Structure="CORTEX_LEFT"
  elif [ $Hemisphere = "R" ] ; then 
    hemisphere="r"
    Structure="CORTEX_RIGHT"
  fi
  
  #native Mesh Processing
  #Convert and volumetrically register white and pial surfaces makign linear and nonlinear copies, add each to the appropriate spec file
  Types="ANATOMICAL@GRAY_WHITE ANATOMICAL@PIAL"
  i=1
  for Surface in white pial ; do
    Type=`echo "$Types" | cut -d " " -f $i`
    Secondary=`echo "$Type" | cut -d "@" -f 2`
    Type=`echo "$Type" | cut -d "@" -f 1`
    if [ ! $Secondary = $Type ] ; then
      Secondary=`echo " -surface-secondary-type ""$Secondary"`
    else
      Secondary=""
    fi
    mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$WBFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    wb_command -set-structure "$WBFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type $Type$Secondary
    wb_command -surface-apply-affine "$WBFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$FreeSurferFolder"/mri/c_ras.mat "$WBFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    wb_command -add-to-spec-file "$WBFolder"/"$Subject".native.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    i=$(($i+1))
  done
  
  #Create midthickness by averaging white and pial surfaces and use it to make inflated surfacess
  wb_command -surface-average "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii -surf "$WBFolder"/"$Subject"."$Hemisphere".white.native.surf.gii -surf "$WBFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii
    wb_command -set-structure "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii ${Structure} -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS
    wb_command -add-to-spec-file "$WBFolder"/"$Subject".native.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii
    wb_command -surface-generate-inflated "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii -iterations-scale 2.5
    wb_command -add-to-spec-file "$WBFolder"/"$Subject".native.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".inflated.native.surf.gii
    wb_command -add-to-spec-file "$WBFolder"/"$Subject".native.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".very_inflated.native.surf.gii
  
  #Convert original and registered spherical surfaces and add them to the nonlinear spec file
  for Surface in sphere.reg sphere ; do
    mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$WBFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
    wb_command -set-structure "$WBFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type SPHERICAL
  done
  wb_command -add-to-spec-file "$WBFolder"/"$Subject".native.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".sphere.native.surf.gii
   
  #Add more files to the spec file and convert other FreeSurfer surface data to metric/GIFTI including sulc, curv, and thickness.
  for Map in sulc@sulc@Sulc thickness@thickness@Thickness curv@curvature@Curvature ; do
    fsname=`echo $Map | cut -d "@" -f 1`
    wbname=`echo $Map | cut -d "@" -f 2`
    mapname=`echo $Map | cut -d "@" -f 3`
    mris_convert -c "$FreeSurferFolder"/surf/"$hemisphere"h."$fsname" "$FreeSurferFolder"/surf/"$hemisphere"h.white "$WBFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
    wb_command -set-structure "$WBFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii ${Structure}
    wb_command -set-map-names "$WBFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -map 1 "$Subject"_"$Hemisphere"_"$mapname"
    wb_command -metric-palette "$WBFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
  done

  #Thickness specific operations
  wb_command -metric-math "abs(thickness)" "$WBFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii -var thickness "$WBFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
  wb_command -metric-palette "$WBFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
  wb_command -metric-math "thickness > 0" "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -var thickness "$WBFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
  wb_command -metric-fill-holes "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
  wb_command -metric-remove-islands "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
  wb_command -set-map-names "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -map 1 "$Subject"_"$Hemisphere"_ROI
  wb_command -metric-dilate "$WBFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$WBFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii -nearest
  wb_command -metric-dilate "$WBFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$WBFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii -nearest

  #Label operations
  for Map in aparc aparc.a2009s BA ; do
    if [ -e "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot ] ; then
      mris_convert --annot "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot "$FreeSurferFolder"/surf/"$hemisphere"h.white "$WBFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii
      wb_command -set-structure "$WBFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii $Structure
      wb_command -set-map-names "$WBFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii -map 1 "$Subject"_"$Hemisphere"_"$Map"
      wb_command -gifti-label-add-prefix "$WBFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii "${Hemisphere}_" "$WBFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii
    fi
  done
  #End main native mesh processing

  #Copy Atlas Files
  cp "$SurfaceAtlasDIR"/fs_"$Hemisphere"/fsaverage."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$WBFolder"/fsaverage/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii
  cp "$SurfaceAtlasDIR"/fs_"$Hemisphere"/fs_"$Hemisphere"-to-fs_LR_fsaverage."$Hemisphere"_LR.spherical_std."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$WBFolder"/fsaverage/"$Subject"."$Hemisphere".def_sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii
  cp "$SurfaceAtlasDIR"/fsaverage."$Hemisphere"_LR.spherical_std."$HighResMesh"k_fs_LR.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii
  wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii
  cp "$SurfaceAtlasDIR"/"$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii
  cp "$SurfaceAtlasDIR"/"$Hemisphere".refsulc."$HighResMesh"k_fs_LR.shape.gii "$WBFolder"/${Subject}.${Hemisphere}.refsulc."$HighResMesh"k_fs_LR.shape.gii
  if [ -e "$SurfaceAtlasDIR"/colin.cerebral."$Hemisphere".flat."$HighResMesh"k_fs_LR.surf.gii ] ; then
    cp "$SurfaceAtlasDIR"/colin.cerebral."$Hemisphere".flat."$HighResMesh"k_fs_LR.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".flat."$HighResMesh"k_fs_LR.surf.gii
    wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".flat."$HighResMesh"k_fs_LR.surf.gii
  fi
  
  #Concatenate FS registration to FS --> FS_LR registration  
  wb_command -surface-sphere-project-unproject "$WBFolder"/"$Subject"."$Hemisphere".sphere.reg.native.surf.gii "$WBFolder"/fsaverage/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$WBFolder"/fsaverage/"$Subject"."$Hemisphere".def_sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$WBFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.surf.gii

  #Make FreeSurfer Registration Areal Distortion Maps
  wb_command -surface-vertex-areas "$WBFolder"/"$Subject"."$Hemisphere".sphere.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".sphere.native.shape.gii
  wb_command -surface-vertex-areas "$WBFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.shape.gii
  wb_command -metric-math "ln(spherereg / sphere) / ln(2)" "$WBFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii -var sphere "$WBFolder"/"$Subject"."$Hemisphere".sphere.native.shape.gii -var spherereg "$WBFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.shape.gii
  rm "$WBFolder"/"$Subject"."$Hemisphere".sphere.native.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.shape.gii
  wb_command -set-map-names "$WBFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii -map 1 "$Subject"_"$Hemisphere"_Areal_Distortion_FS
  wb_command -metric-palette "$WBFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii MODE_AUTO_SCALE -palette-name ROY-BIG-BL -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_OUTSIDE -1 1

   RegSphere="${WBFolder}/${Subject}.${Hemisphere}.sphere.reg.reg_LR.native.surf.gii"

  ##Ensure no zeros in atlas medial wall ROI
  wb_command -metric-resample "$WBFolder"/"$Subject"."$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ${RegSphere} BARYCENTRIC "$WBFolder"/"$Subject"."$Hemisphere".atlasroi.native.shape.gii -largest
  wb_command -metric-math "(atlas + individual) > 0" "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -var atlas "$WBFolder"/"$Subject"."$Hemisphere".atlasroi.native.shape.gii -var individual "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
  wb_command -metric-mask "$WBFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
  wb_command -metric-mask "$WBFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii


  #Populate Highres fs_LR spec file.  Deform surfaces and other data according to native to folding-based registration selected above.  Regenerate inflated surfaces.
  for Surface in white midthickness pial ; do
    wb_command -surface-resample "$WBFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii BARYCENTRIC "$WBFolder"/"$Subject"."$Hemisphere"."$Surface"."$HighResMesh"k_fs_LR.surf.gii
    wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere"."$Surface"."$HighResMesh"k_fs_LR.surf.gii
  done
  wb_command -surface-generate-inflated "$WBFolder"/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".inflated."$HighResMesh"k_fs_LR.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".very_inflated."$HighResMesh"k_fs_LR.surf.gii -iterations-scale 2.5
  wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".inflated."$HighResMesh"k_fs_LR.surf.gii
  wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".very_inflated."$HighResMesh"k_fs_LR.surf.gii
  
  for Map in thickness curvature ; do
    wb_command -metric-resample "$WBFolder"/"$Subject"."$Hemisphere"."$Map".native.shape.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$WBFolder"/"$Subject"."$Hemisphere"."$Map"."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii -current-roi "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
    wb_command -metric-mask "$WBFolder"/"$Subject"."$Hemisphere"."$Map"."$HighResMesh"k_fs_LR.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii "$WBFolder"/"$Subject"."$Hemisphere"."$Map"."$HighResMesh"k_fs_LR.shape.gii
  done  
  wb_command -metric-resample "$WBFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$WBFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii
  wb_command -metric-resample "$WBFolder"/"$Subject"."$Hemisphere".sulc.native.shape.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$WBFolder"/"$Subject"."$Hemisphere".sulc."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii

  for Map in aparc aparc.a2009s BA ; do
    if [ -e "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot ] ; then
      wb_command -label-resample "$WBFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii BARYCENTRIC "$WBFolder"/"$Subject"."$Hemisphere"."$Map"."$HighResMesh"k_fs_LR.label.gii -largest
    fi
  done

  for LowResMesh in ${LowResMeshes} ; do
    #Copy Atlas Files
    cp "$SurfaceAtlasDIR"/"$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii
    wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii
    cp "$SurfaceAtlasDIR"/"$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii
    if [ -e "$SurfaceAtlasDIR"/colin.cerebral."$Hemisphere".flat."$LowResMesh"k_fs_LR.surf.gii ] ; then
      cp "$SurfaceAtlasDIR"/colin.cerebral."$Hemisphere".flat."$LowResMesh"k_fs_LR.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".flat."$LowResMesh"k_fs_LR.surf.gii
      wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".flat."$LowResMesh"k_fs_LR.surf.gii
    fi

    #Create downsampled fs_LR spec files.   
    for Surface in white midthickness pial ; do
      wb_command -surface-resample "$WBFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii BARYCENTRIC "$WBFolder"/"$Subject"."$Hemisphere"."$Surface"."$LowResMesh"k_fs_LR.surf.gii
      wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere"."$Surface"."$LowResMesh"k_fs_LR.surf.gii
    done

    wb_command -surface-generate-inflated "$WBFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".inflated."$LowResMesh"k_fs_LR.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".very_inflated."$LowResMesh"k_fs_LR.surf.gii -iterations-scale 0.75
    wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".inflated."$LowResMesh"k_fs_LR.surf.gii
    wb_command -add-to-spec-file "$WBFolder"/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$WBFolder"/"$Subject"."$Hemisphere".very_inflated."$LowResMesh"k_fs_LR.surf.gii
  
    for Map in sulc thickness curvature ; do
      wb_command -metric-resample "$WBFolder"/"$Subject"."$Hemisphere"."$Map".native.shape.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$WBFolder"/"$Subject"."$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii -current-roi "$WBFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
      wb_command -metric-mask "$WBFolder"/"$Subject"."$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.shape.gii "$WBFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii "$WBFolder"/"$Subject"."$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.shape.gii
    done  
    wb_command -metric-resample "$WBFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$WBFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii
    wb_command -metric-resample "$WBFolder"/"$Subject"."$Hemisphere".sulc.native.shape.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$WBFolder"/"$Subject"."$Hemisphere".sulc."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$WBFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$WBFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii

    for Map in aparc aparc.a2009s BA ; do
      if [ -e "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot ] ; then
        wb_command -label-resample "$WBFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii ${RegSphere} "$WBFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii BARYCENTRIC "$WBFolder"/"$Subject"."$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.label.gii -largest
      fi
    done

  done
done

#STRINGII=""
#for LowResMesh in ${LowResMeshes} ; do
#  STRINGII=`echo "${STRINGII}${AtlasSpaceFolder}/fsaverage${LowResMesh}k@${LowResMesh}k_fs_LR@atlasroi "`
#done
#
##Create CIFTI Files
#for STRING in "$AtlasSpaceFolder"/"$NativeFolder"@native@roi "$AtlasSpaceFolder"@"$HighResMesh"k_fs_LR@atlasroi ${STRINGII} ; do
#  Folder=`echo $STRING | cut -d "@" -f 1`
#  Mesh=`echo $STRING | cut -d "@" -f 2`
#  ROI=`echo $STRING | cut -d "@" -f 3`
#  
#  wb_command -cifti-create-dense-scalar "$Folder"/"$Subject".sulc."$Mesh".dscalar.nii -left-metric "$Folder"/"$Subject".L.sulc."$Mesh".shape.gii -right-metric "$Folder"/"$Subject".R.sulc."$Mesh".shape.gii
#  wb_command -set-map-names "$Folder"/"$Subject".sulc."$Mesh".dscalar.nii -map 1 "${Subject}_Sulc"
#  wb_command -cifti-palette "$Folder"/"$Subject".sulc."$Mesh".dscalar.nii MODE_AUTO_SCALE_PERCENTAGE "$Folder"/"$Subject".sulc."$Mesh".dscalar.nii -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
#      
#  wb_command -cifti-create-dense-scalar "$Folder"/"$Subject".curvature."$Mesh".dscalar.nii -left-metric "$Folder"/"$Subject".L.curvature."$Mesh".shape.gii -roi-left "$Folder"/"$Subject".L."$ROI"."$Mesh".shape.gii -right-metric "$Folder"/"$Subject".R.curvature."$Mesh".shape.gii -roi-right "$Folder"/"$Subject".R."$ROI"."$Mesh".shape.gii
#  wb_command -set-map-names "$Folder"/"$Subject".curvature."$Mesh".dscalar.nii -map 1 "${Subject}_Curvature"
#  wb_command -cifti-palette "$Folder"/"$Subject".curvature."$Mesh".dscalar.nii MODE_AUTO_SCALE_PERCENTAGE "$Folder"/"$Subject".curvature."$Mesh".dscalar.nii -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
#
#  wb_command -cifti-create-dense-scalar "$Folder"/"$Subject".thickness."$Mesh".dscalar.nii -left-metric "$Folder"/"$Subject".L.thickness."$Mesh".shape.gii -roi-left "$Folder"/"$Subject".L."$ROI"."$Mesh".shape.gii -right-metric "$Folder"/"$Subject".R.thickness."$Mesh".shape.gii -roi-right "$Folder"/"$Subject".R."$ROI"."$Mesh".shape.gii
#  wb_command -set-map-names "$Folder"/"$Subject".thickness."$Mesh".dscalar.nii -map 1 "${Subject}_Thickness"
#  wb_command -cifti-palette "$Folder"/"$Subject".thickness."$Mesh".dscalar.nii MODE_AUTO_SCALE_PERCENTAGE "$Folder"/"$Subject".thickness."$Mesh".dscalar.nii -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
# 
#  wb_command -cifti-create-dense-scalar "$Folder"/"$Subject".ArealDistortion_FS."$Mesh".dscalar.nii -left-metric "$Folder"/"$Subject".L.ArealDistortion_FS."$Mesh".shape.gii -right-metric "$Folder"/"$Subject".R.ArealDistortion_FS."$Mesh".shape.gii
#  wb_command -set-map-names "$Folder"/"$Subject".ArealDistortion_FS."$Mesh".dscalar.nii -map 1 "${Subject}_ArealDistortion_FS"
#  wb_command -cifti-palette "$Folder"/"$Subject".ArealDistortion_FS."$Mesh".dscalar.nii MODE_USER_SCALE "$Folder"/"$Subject".ArealDistortion_FS."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false
#
#  if [ ${RegName} = "MSMSulc" ] ; then
#    wb_command -cifti-create-dense-scalar "$Folder"/"$Subject".ArealDistortion_MSMSulc."$Mesh".dscalar.nii -left-metric "$Folder"/"$Subject".L.ArealDistortion_MSMSulc."$Mesh".shape.gii -right-metric "$Folder"/"$Subject".R.ArealDistortion_MSMSulc."$Mesh".shape.gii
#    wb_command -set-map-names "$Folder"/"$Subject".ArealDistortion_MSMSulc."$Mesh".dscalar.nii -map 1 "${Subject}_ArealDistortion_MSMSulc"
#    wb_command -cifti-palette "$Folder"/"$Subject".ArealDistortion_MSMSulc."$Mesh".dscalar.nii MODE_USER_SCALE "$Folder"/"$Subject".ArealDistortion_MSMSulc."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false
#  fi
#  
#  for Map in aparc aparc.a2009s BA ; do 
#    if [ -e "$Folder"/"$Subject".L.${Map}."$Mesh".label.gii ] ; then
#      wb_command -cifti-create-label "$Folder"/"$Subject".${Map}."$Mesh".dlabel.nii -left-label "$Folder"/"$Subject".L.${Map}."$Mesh".label.gii -roi-left "$Folder"/"$Subject".L."$ROI"."$Mesh".shape.gii -right-label "$Folder"/"$Subject".R.${Map}."$Mesh".label.gii -roi-right "$Folder"/"$Subject".R."$ROI"."$Mesh".shape.gii
#      wb_command -set-map-names "$Folder"/"$Subject".${Map}."$Mesh".dlabel.nii -map 1 "$Subject"_${Map}
#    fi
#  done 
#done
#
#STRINGII=""
#for LowResMesh in ${LowResMeshes} ; do
#  STRINGII=`echo "${STRINGII}${AtlasSpaceFolder}/fsaverage${LowResMesh}k@${AtlasSpaceFolder}/fsaverage${LowResMesh}k@${LowResMesh}k_fs_LR ${T1wFolder}/fsaverage${LowResMesh}k@${AtlasSpaceFolder}/fsaveragfsaveragfsaveragfsaverage${LowResMesh}k@${LowResMesh}k_fs_LR "`
#done
#
##Add CIFTI Maps to Spec Files
#for STRING in "$T1wFolder"/"$NativeFolder"@"$AtlasSpaceFolder"/"$NativeFolder"@native "$AtlasSpaceFolder"/"$NativeFolder"@"$AtlasSpaceFolder"/"$NativeFolder"@native "$AtlasSpaceFolder"@"$AtlasSpaceFolder"@"$HighResMesh"k_fs_LR ${STRINGII} ; do
#  FolderI=`echo $STRING | cut -d "@" -f 1`
#  FolderII=`echo $STRING | cut -d "@" -f 2`
#  Mesh=`echo $STRING | cut -d "@" -f 3`
#  for STRINGII in sulc@dscalar thickness@dscalar curvature@dscalar aparc@dlabel aparc.a2009s@dlabel BA@dlabel ; do
#    Map=`echo $STRINGII | cut -d "@" -f 1`
#    Ext=`echo $STRINGII | cut -d "@" -f 2`
#    if [ -e "$FolderII"/"$Subject"."$Map"."$Mesh"."$Ext".nii ] ; then
#      wb_command -add-to-spec-file "$FolderI"/"$Subject"."$Mesh".wb.spec INVALID "$FolderII"/"$Subject"."$Map"."$Mesh"."$Ext".nii
#    fi
#  done
#done

echo -e "\n END: FS2CaretConvertRegisterNonlinear"



