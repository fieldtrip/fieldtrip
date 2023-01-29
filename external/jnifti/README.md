# JNIfTI Toolbox - Fast and portable NIfTI-1/2 reader and NIfTI-to-JNIfTI converter

* Copyright (C) 2019  Qianqian Fang <q.fang at neu.edu>
* License: GNU General Public License version 3 (GPL v3) or Apache License 2.0, see License*.txt
* Version: 0.5 (Ascendence)
* URL: http://github.com/fangq/jnifti

## Overview

This is a fully functional NIfTI-1/2 reader/writer that supports both 
MATLAB and GNU Octave, and is capable of reading/writing both non-compressed 
and compressed NIfTI files (.nii, .nii.gz) as well as two-part Analyze7.5/NIfTI
files (.hdr/.img and .hdr.gz/.img.gz).

More importantly, this is a toolbox that converts NIfTI data to its JSON-based
replacement, JNIfTI (.jnii for text-based and .bnii for binary-based), defined
by the JNIfTI specification (http://github.com/fangq/jnifti). JNIfTI is a 
much more flexible, human-readable and extensible file format compared to the
more rigid and opaque NIfTI format, making the data much easier to manipulate
and share.

## Installation

The JNIfTI toolbox includes a stand-alone NIfTI-1/2 parser that works on both
MATLAB and GNU Octave without needing additional components. To just reading and
writing the un-compressed NIfTI and Analyze7.5 files (.nii, .hdr/.img), one 
only needs to run `addpath('/path/to/jnifti')`. For MATLAB, JNIfTI toolbox
utilizes `memmapfile`-based disk-reading, making it very fast. For Octave, 
`memmapfile` is currently not implemented, so, a full reading is required.

The JNIfTI toolbox is also capable of reading/writing gzip-compressed NIfTI and 
Analyze7.5 files (.nii.gz, .hdr.gz, .img.gz). This feature is supported in MATLAB
directly without needing another toolbox (MATLAB must be in the JVM-enabled mode).

To process gzip-compressed NIfTI/Analyze files in Octave and MATLAB with `-nojvm`,
one need to install the open-source JSONLab and ZMat toolboxes, both supporting
MATLAB and Octave. They can be downloaded at

* JSONLab: http://github.com/fangq/jsonlab
* ZMat: http://github.com/fangq/zmat

To save NIfTI-1/2 data as JNIfTI files, one needs to install JSONLab. The JNIfTI
data format supports internal compression (as oppose to external compression such
as \*.gz files). To create or read compressed JNIfTI files in Octave, one must 
install the ZMat toolbox, as listed above.

## Usage

### `nii2jnii` - To convert a NIfTI-1/2 file to a JNIfTI file or data structure
Example:
```
  nii=nii2jnii('test.nii', 'nii');         % read a .nii file as a nii structure
  nii=nii2jnii('test.nii.gz');             % read a .nii.gz file as a jnii structure
  nii2jnii('test.nii.gz', 'newdata.jnii') ;% read a .nii.gz file and convert to a text-JNIfTI file  
  nii2jnii('test.nii.gz', 'newdata.bnii','compression','zlib'); % read a .nii.gz file and convert to a binary-JNIfTI file with compression
```
### `loadnifti` - To read a NIfTI-1/2 (.nii or .nii.gz) file (alias to `nii2jnii`)
Example:
```
  nii=loadnifti('test.nii.gz');             % read a .nii.gz file as a jnii structure
  nii=loadnifti('test.nii', 'nii');         % read a .nii file as a nii structure
```
### `savenifti` - To write an image as NIfTI-1/2 (.nii or .nii.gz) file
Example:
```
  savenifti(img,'test.nii.gz');         % save an array img to a compressed nifti file
  savenifti(img, 'test.nii', 'nifti2'); % save an array img to a nifti-2 file file
  savenifti(img, 'test.nii', header);   % save an array img with an existing header
```
### `loadjnifti` - To read a JNIfTI (.jnii or .bnii) file
Example:
```
  jnii=nii2jnii('test.nii.gz');
  savejnifti(jnii, 'magic10.bnii','Compression','gzip');
  newjnii=loadjnifti('magic10.bnii');
```
### `savejnifti` - To write a JNIfTI structure into a file (.jnii or .bnii)
Example:
```
  jnii=jnifticreate(uint8(magic(10)),'Name','10x10 magic matrix');
  savejnifti(jnii, 'magic10.jnii');
  savejnifti(jnii, 'magic10_debug.bnii','Compression','gzip');
```
