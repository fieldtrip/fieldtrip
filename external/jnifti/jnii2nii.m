function nii=jnii2nii(jnii, varargin)
%
%    nii=jnii2nii(jnii)
%       or
%    nii=jnii2nii(jniifile)
%    jnii2nii(jniifile, niifile)
%
%    Covert a JNIfTI file or data structure to a NIfTI-1/2 structure or file
%
%    This function is compatible with both MATLAB and GNU Octave. 
%    It accepts .jnii and .bnii input files
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        jnii: a JNIfTI data structure (a struct with NIFTIHeader and NIFTIData fields);
%              if jnii is a string, it represents a JNIfTI file (.jnii/.bnii)
%        niifile: if the 2nd parameter is given as a file name, the converted nifti data
%              will be save as a nii file with filename specified by niifile.
%              if the filename in niifile contains .gz, the file will be compressed using
%              the zmat toolbox.
%
%    output:
%        nii: is the converted nifti-1/2 data structure, it contains the below subfields
%          nii.img: the data volume read from the nii file
%          nii.hdr: extended raw file header, a structure that is byte-wise compatible with a 
%                   nifti-1 - in this case, typecast(nii.hdr,'uint8') must be 348+4=352 bytes, 
%                       including the raw nifti-1 hdr header (348 bytes) plus the 4-byte 
%                       extension flags), or
%                   nifti-2 - in this case, typecast(nii.hdr,'uint8') must be 540+4=544 bytes, 
%                       including the raw nifti-2 hdr header (540 bytes) plus the 4-byte 
%                       extension flags)
%              if one run nii.hdr.extension=[]; the resulting struct is 348/540-byte in length
%              nii.hdr key subfileds include
%
%              sizeof_hdr: must be 348 (for NIFTI-1) or 540 (for NIFTI-2)
%              dim: short array, dim(2: dim(1)+1) defines the array size
%              datatype: the type of data stored in each voxel
%              bitpix: total bits per voxel
%              magic: must be 'ni1\0' or 'n+1\0' for NIFTI-1 data, and 'ni2\0' or 'n+2\0' for NIFTI-2 data
%
%              For the detailed nii header, please see 
%              https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
%
%    dependency:
%
%          To load a JNIfTI file with compression or niifile ends with (.nii.gz/.hdr.gz/.img.gz),
%          one must install the ZMat Toolbox (http://github.com/fangq/zmat) and
%          JSONLab Toolbox (http://github.com/fangq/jsonlab);
%
%    this file is part of JNIfTI specification: https://github.com/fangq/jnifti
%
%    License: Apache 2.0, see https://github.com/fangq/jnifti for details
%

if(nargin<=0)
    help jnii2nii
    return;
end

if(~isstruct(jnii))
    jnii=loadjnifti(jnii);
end

if(~(isfield(jnii,'NIFTIHeader') && isfield(jnii,'NIFTIData')))
    error('input must be a valid JNIfTI structure (needs both NIFTIHeader and NIFTIData subfields)');
end

niiformat='nifti1';

if((isfield(jnii.NIFTIHeader,'NIIFormat') && ismember(jnii.NIFTIHeader.NIIFormat(1:3),{'ni2','n+2'})) || max(jnii.NIFTIHeader.Dim)>=2^32)
    niiformat='nifti2';
end

nii.hdr=nifticreate(jnii.NIFTIData, niiformat);
nii.img=jnii.NIFTIData;

if(isfield(jnii.NIFTIHeader,'NIIHeaderSize'))
    nii.hdr.sizeof_hdr	=bytematch(jnii.NIFTIHeader, 'NIIHeaderSize', nii.hdr.sizeof_hdr);
end
if(isfield(nii.hdr,'data_type'))
    nii.hdr.data_type	=bytematch(jnii.NIFTIHeader, 'A75DataTypeName', nii.hdr.data_type);
    nii.hdr.db_name	=bytematch(jnii.NIFTIHeader, 'A75DBName', nii.hdr.db_name);
    nii.hdr.extents	=bytematch(jnii.NIFTIHeader, 'A75Extends', nii.hdr.extents);
    nii.hdr.session_error	=bytematch(jnii.NIFTIHeader, 'A75SessionError', nii.hdr.session_error);
    nii.hdr.regular	=bytematch(jnii.NIFTIHeader, 'A75Regular', nii.hdr.regular);
end

dim_info=bitor(uint8(jnii.NIFTIHeader.DimInfo.Freq),bitshift(uint8(jnii.NIFTIHeader.DimInfo.Phase),3));
dim_info=bitor(dim_info,bitshift(uint8(jnii.NIFTIHeader.DimInfo.Slice),6));
nii.hdr.dim_info=cast(dim_info, class(nii.hdr.dim_info));

nii.hdr.dim(1)=cast(length(jnii.NIFTIHeader.Dim),class(nii.hdr.dim));
nii.hdr.dim(2:1+length(jnii.NIFTIHeader.Dim))	=bytematch(jnii.NIFTIHeader, 'Dim', nii.hdr.dim(2:1+length(jnii.NIFTIHeader.Dim)));
nii.hdr.intent_p1	=bytematch(jnii.NIFTIHeader, 'Param1', nii.hdr.intent_p1);
nii.hdr.intent_p2	=bytematch(jnii.NIFTIHeader, 'Param2', nii.hdr.intent_p2);
nii.hdr.intent_p3	=bytematch(jnii.NIFTIHeader, 'Param3', nii.hdr.intent_p3);

if(ischar(jnii.NIFTIHeader.Intent))
    jnii.NIFTIHeader.Intent=niicodemap('intent',jnii.NIFTIHeader.Intent);
end
nii.hdr.intent_code	=bytematch(jnii.NIFTIHeader, 'Intent', nii.hdr.intent_code);
if(ischar(jnii.NIFTIHeader.DataType))
    jnii.NIFTIHeader.DataType=niicodemap('datatype',jnii.NIFTIHeader.DataType);
end

nii.hdr.datatype	=bytematch(jnii.NIFTIHeader, 'DataType', nii.hdr.datatype);
nii.hdr.bitpix	=bytematch(jnii.NIFTIHeader, 'BitDepth', nii.hdr.bitpix);
nii.hdr.slice_start	=bytematch(jnii.NIFTIHeader, 'FirstSliceID', nii.hdr.slice_start);
nii.hdr.pixdim(1)=cast(length(jnii.NIFTIHeader.VoxelSize),class(nii.hdr.pixdim));
nii.hdr.pixdim(2:2+nii.hdr.dim(1)-1)	=bytematch(jnii.NIFTIHeader, 'VoxelSize', nii.hdr.pixdim(2:2+nii.hdr.dim(1)-1));
nii.hdr.vox_offset	=bytematch(jnii.NIFTIHeader, 'NIIByteOffset', nii.hdr.vox_offset);
nii.hdr.scl_slope	=bytematch(jnii.NIFTIHeader, 'ScaleSlope', nii.hdr.scl_slope);
nii.hdr.scl_inter	=bytematch(jnii.NIFTIHeader, 'ScaleOffset', nii.hdr.scl_inter);
nii.hdr.slice_end	=bytematch(jnii.NIFTIHeader, 'LastSliceID', nii.hdr.slice_end);

if(ischar(jnii.NIFTIHeader.SliceType))
    jnii.NIFTIHeader.SliceType=niicodemap('slicetype',jnii.NIFTIHeader.SliceType);
end
nii.hdr.slice_code	=bytematch(jnii.NIFTIHeader, 'SliceType', nii.hdr.slice_code);
if(ischar(jnii.NIFTIHeader.Unit.L))
    jnii.NIFTIHeader.Unit.L=niicodemap('unit',jnii.NIFTIHeader.Unit.L);
end
if(ischar(jnii.NIFTIHeader.Unit.T))
    jnii.NIFTIHeader.Unit.T=niicodemap('unit',jnii.NIFTIHeader.Unit.T);
end

xyzt_units=bitor(uint8(jnii.NIFTIHeader.Unit.L),uint8(jnii.NIFTIHeader.Unit.T));
nii.hdr.xyzt_units=cast(xyzt_units, class(nii.hdr.xyzt_units));

nii.hdr.cal_max	=bytematch(jnii.NIFTIHeader, 'MaxIntensity', nii.hdr.cal_max);
nii.hdr.cal_min	=bytematch(jnii.NIFTIHeader, 'MinIntensity', nii.hdr.cal_min);
nii.hdr.slice_duration	=bytematch(jnii.NIFTIHeader, 'SliceTime', nii.hdr.slice_duration);
nii.hdr.toffset	=bytematch(jnii.NIFTIHeader, 'TimeOffset', nii.hdr.toffset);
if(isfield(nii.hdr,'glmax'))
    nii.hdr.glmax	=bytematch(jnii.NIFTIHeader, 'A75GlobalMax', nii.hdr.glmax);
    nii.hdr.glmin	=bytematch(jnii.NIFTIHeader, 'A75GlobalMin', nii.hdr.glmin);
end

nii.hdr.descrip	=bytematch(jnii.NIFTIHeader, 'Description', nii.hdr.descrip);
nii.hdr.aux_file	=bytematch(jnii.NIFTIHeader, 'AuxFile', nii.hdr.aux_file);

nii.hdr.qform_code	=bytematch(jnii.NIFTIHeader, 'QForm', nii.hdr.qform_code);
nii.hdr.sform_code	=bytematch(jnii.NIFTIHeader, 'SForm', nii.hdr.sform_code);
nii.hdr.quatern_b	=bytematch(jnii.NIFTIHeader, 'Quatern.b', nii.hdr.quatern_b);
nii.hdr.quatern_c	=bytematch(jnii.NIFTIHeader, 'Quatern.c', nii.hdr.quatern_c);
nii.hdr.quatern_d	=bytematch(jnii.NIFTIHeader, 'Quatern.d', nii.hdr.quatern_d);
nii.hdr.qoffset_x	=bytematch(jnii.NIFTIHeader, 'QuaternOffset.x', nii.hdr.qoffset_x);
nii.hdr.qoffset_y	=bytematch(jnii.NIFTIHeader, 'QuaternOffset.y', nii.hdr.qoffset_y);
nii.hdr.qoffset_z	=bytematch(jnii.NIFTIHeader, 'QuaternOffset.z', nii.hdr.qoffset_z);
nii.hdr.srow_x	=cast(jnii.NIFTIHeader.Affine(1,:), class(nii.hdr.srow_x));
nii.hdr.srow_y	=cast(jnii.NIFTIHeader.Affine(2,:), class(nii.hdr.srow_y));
nii.hdr.srow_z	=cast(jnii.NIFTIHeader.Affine(3,:), class(nii.hdr.srow_z));

nii.hdr.intent_name	=bytematch(jnii.NIFTIHeader, 'Name', nii.hdr.intent_name);
%nii.hdr.magic	=bytematch(jnii.NIFTIHeader, 'NIIFormat', nii.hdr.magic);

if(isfield(jnii.NIFTIHeader,'NIIExtender'))
    nii.hdr.extension	=bytematch(jnii.NIFTIHeader, 'NIIExtender', nii.hdr.extension);
end
if(isfield(jnii.NIFTIHeader,'NIIQfac_'))
    nii.hdr.pixdim(1)	=bytematch(jnii.NIFTIHeader, 'NIIQfac_', nii.hdr.pixdim(1));
end
if(isfield(jnii.NIFTIHeader,'NIIUnused_'))
    nii.hdr.reserved	=bytematch(jnii.NIFTIHeader, 'NIIUnused_', nii.hdr.reserved);
end

if(isfield(jnii,'NIFTIExtension') && iscell(jnii.NIFTIExtension))
    nii.extension=jnii.NIFTIExtension;
    if(nii.hdr.extension(1)~=length(jnii.NIFTIExtension))
        nii.hdr.extension(1)=length(jnii.NIFTIExtension);
	warning('header extension count does not match the extension data, force update');
    end
end

if(nargin>=2 && ischar(varargin{1}))
    savenifti(nii.img, varargin{1}, nii.hdr);
end

%---------------------------------------------------------------------------

function dat=bytematch(jobj, key, orig)
dtype=class(orig);
dat=orig;
if(isfield(jobj,key))
    dat=cast(jobj.(key),dtype);
else
    dat=cast(0,dtype);
end
if(length(dat)<length(orig))
    dat(length(orig))=cast(0,dtype);
end
dat=dat(1:length(orig));
