function nii=readnifti(filename)
%
%    nii=readnifti(filename)
%
%    Read a Nifti (*.nii) or Analyze 7.5 (*.hdr/*.img) image file
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        fname: the file name to a .nii file, or an Analyze 7.5 file (*.hdr,*.img)
%
%    output:
%        nii.img: the data volume read from the nii file
%        nii.datatype: the data type of the voxel, in matlab data type string
%        nii.datalen: data count per voxel - for example RGB data has 3x
%                    uint8 per voxel, so datatype='uint8', datalen=3
%        nii.voxelbyte: total number of bytes per voxel: for RGB data,
%                    voxelbyte=3; also voxelbyte=header.bitpix/8
%        nii.hdr: file header info, a structure has the full nii header
%                    key subfileds include
%
%              sizeof_hdr: must be 348 if the input is nifti
%              dim: short array, dim(2: dim(1)+1) defines the array size
%              datatype: the type of data stored in each voxel
%              bitpix: total bits per voxel
%              magic: must be 'ni1\0' or 'n+1\0'
%
%              For the detailed nii header, please see
%              https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
%
%   this file was ported from mcxloadnii.m from the MCX Project (http://mcx.space)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

hdrfile=filename;
isnii=-1;
if(regexp(filename,'(\.[Hh][Dd][Rr]$|\.[Ii][Mm][Gg]$)'))
    isnii=0;
elseif(regexp(filename,'\.[Nn][Ii][Ii]$'))
    isnii=1;
end

if(isnii<0)
    error('file is not a Nifti or Analyze 7.5 data file');
end

if(regexp(filename,'\.[Ii][Mm][Gg]$'))
    hdrfile=regexprep(filename,'\.[Ii][Mm][Gg]$','.hdr');
end

header = memmapfile(hdrfile,                ...
      'Offset', 0,                           ...
      'Writable', false,                     ...
      'Format', {                            ...
     'int32'   [1 1]  'sizeof_hdr'    ; %!< MUST be 348 	      %  % int sizeof_hdr;       %  ...
     'int8'    [1 10] 'data_type'     ; %!< ++UNUSED++  	      %  % char data_type[10];   %  ...
     'int8'    [1 18] 'db_name'       ; %!< ++UNUSED++  	      %  % char db_name[18];     %  ...
     'int32'   [1 1]  'extents'       ; %!< ++UNUSED++  	      %  % int extents;	         %  ...
     'int16'   [1 1]  'session_error' ; %!< ++UNUSED++  	      %  % short session_error;  %  ...
     'int8'    [1 1]  'regular'       ; %!< ++UNUSED++  	      %  % char regular;	     %  ...
     'int8'    [1 1]  'dim_info'      ; %!< MRI slice ordering.   %  % char hkey_un0;	     %  ...
     'uint16'  [1 8]  'dim'	          ; %!< Data array dimensions.%  % short dim[8];	     %  ...
     'single'  [1 1]  'intent_p1'     ; %!< 1st intent parameter. %  % short unused8/9;      %  ...
     'single'  [1 1]  'intent_p2'     ; %!< 2nd intent parameter. %  % short unused10/11;    %  ...
     'single'  [1 1]  'intent_p3'     ; %!< 3rd intent parameter. %  % short unused12/13;    %  ...
     'int16'   [1 1]  'intent_code'   ; %!< NIFTI_INTENT_* code.  %  % short unused14;       %  ...
     'int16'   [1 1]  'datatype'      ; %!< Defines data type!    %  % short datatype;       %  ...
     'int16'   [1 1]  'bitpix'        ; %!< Number bits/voxel.    %  % short bitpix;	     %  ...
     'int16'   [1 1]  'slice_start'   ; %!< First slice index.    %  % short dim_un0;	     %  ...
     'single'  [1 8]  'pixdim'        ; %!< Grid spacings.	      %  % float pixdim[8];      %  ...
     'single'  [1 1]  'vox_offset'    ; %!< Offset into .nii file %  % float vox_offset;     %  ...
     'single'  [1 1]  'scl_slope'     ; %!< Data scaling: slope.  %  % float funused1;       %  ...
     'single'  [1 1]  'scl_inter'     ; %!< Data scaling: offset. %  % float funused2;       %  ...
     'int16'   [1 1]  'slice_end'     ; %!< Last slice index.	  %  % float funused3;       %  ...
     'int8'    [1 1]  'slice_code'    ; %!< Slice timing order.   %				                ...
     'int8'    [1 1]  'xyzt_units'    ; %!< Units of pixdim[1..4] %				                ...
     'single'  [1 1]  'cal_max'       ; %!< Max display intensity %  % float cal_max;	     %  ...
     'single'  [1 1]  'cal_min'       ; %!< Min display intensity %  % float cal_min;	     %  ...
     'single'  [1 1]  'slice_duration'; %!< Time for 1 slice.	  %  % float compressed;     %  ...
     'single'  [1 1]  'toffset'       ; %!< Time axis shift.	  %  % float verified;       %  ...
     'int32'   [1 1]  'glmax'	      ; %!< ++UNUSED++  	      %  % int glmax;	         %  ...
     'int32'   [1 1]  'glmin'	      ; %!< ++UNUSED++  	      %  % int glmin;	         %  ...
     'int8'    [1 80] 'descrip'       ; %!< any text you like.    %  % char descrip[80];     %  ...
     'int8'    [1 24] 'aux_file'      ; %!< auxiliary filename.   %  % char aux_file[24];    %  ...
     'int16'   [1 1]  'qform_code'    ; %!< NIFTI_XFORM_* code.   %  %-- all ANALYZE 7.5 --- %  ...
     'int16'   [1 1]  'sform_code'    ; %!< NIFTI_XFORM_* code.   %  %below here are replaced%  ...
     'single'  [1 1]  'quatern_b'     ; %!< Quaternion b param.   %				...
     'single'  [1 1]  'quatern_c'     ; %!< Quaternion c param.   %				...
     'single'  [1 1]  'quatern_d'     ; %!< Quaternion d param.   %				...
     'single'  [1 1]  'qoffset_x'     ; %!< Quaternion x shift.   %				...
     'single'  [1 1]  'qoffset_y'     ; %!< Quaternion y shift.   %				...
     'single'  [1 1]  'qoffset_z'     ; %!< Quaternion z shift.   %				...
     'single'  [1 4]  'srow_x'        ; %!< 1st row affine transform.	%			...
     'single'  [1 4]  'srow_y'        ; %!< 2nd row affine transform.	%			...
     'single'  [1 4]  'srow_z'        ; %!< 3rd row affine transform.	%			...
     'int8'    [1 16] 'intent_name'   ; %!< 'name' or meaning of data.  %			...
     'int8'    [1 4]  'magic'	      ; %!< MUST be "ni1\0" or "n+1\0". %			...
     'int8'    [1 4]  'extension'	    %!< header extension	  %				...
         });

nii.hdr=header.Data(1);
     
type2byte=[
        0  0  % unknown                      %
        1  0  % binary (1 bit/voxel)         %
        2  1  % unsigned char (8 bits/voxel) %
        4  2  % signed short (16 bits/voxel) %
        8  4  % signed int (32 bits/voxel)   %
       16  4  % float (32 bits/voxel)        %
       32  8  % complex (64 bits/voxel)      %
       64  8  % double (64 bits/voxel)       %
      128  3  % RGB triple (24 bits/voxel)   %
      255  0  % not very useful (?)          %
      256  1  % signed char (8 bits)         %
      512  2  % unsigned short (16 bits)     %
      768  4  % unsigned int (32 bits)       %
     1024  8  % long long (64 bits)          %
     1280  8  % unsigned long long (64 bits) %
     1536 16  % long double (128 bits)       %
     1792 16  % double pair (128 bits)       %
     2048 32  % long double pair (256 bits)  %
     2304  4  % 4 byte RGBA (32 bits/voxel)  %
];

type2str={
    'uint8'    0   % unknown                       %
    'uint8'    0   % binary (1 bit/voxel)          %
    'uint8'    1   % unsigned char (8 bits/voxel)  %
    'uint16'   1   % signed short (16 bits/voxel)  %
    'int32'    1   % signed int (32 bits/voxel)    %
    'float32'  1   % float (32 bits/voxel)         %
    'float32'  2   % complex (64 bits/voxel)       %
    'float64'  1   % double (64 bits/voxel)        %
    'uint8'    3   % RGB triple (24 bits/voxel)    %
    'uint8'    0   % not very useful (?)           %
    'int8'     1   % signed char (8 bits)          %
    'uint16'   1   % unsigned short (16 bits)      %
    'uint32'   1   % unsigned int (32 bits)        %
    'long'     1   % long long (64 bits)           %
    'ulong'    1   % unsigned long long (64 bits)  %
    'uint8'    16  % long double (128 bits)        %
    'uint8'    16  % double pair (128 bits)        %
    'uint8'    32  % long double pair (256 bits)   %
    'uint8'    4   % 4 byte RGBA (32 bits/voxel)   %
};

typeidx=find(type2byte(:,1)==nii.hdr.datatype);

nii.datatype=type2str{typeidx,1};
nii.datalen=type2str{typeidx,2};
nii.voxelbyte=type2byte(typeidx,2);

if(type2byte(typeidx,2)==0)
    nii.img=[];
    return;
end

if(type2str{typeidx,2}>1)
    nii.hdr.dim=[nii.hdr.dim(1)+1 uint16(nii.datalen) nii.hdr.dim(2:end)]; 
end

if(regexp(filename,'\.[Hh][Dd][Rr]$'))
    filename=regexprep(filename,'\.[Hh][Dd][Rr]$','.img');
end

fid=fopen(filename,'rb');
if(isnii)
   fseek(fid,nii.hdr.vox_offset,'bof');
end
nii.img=fread(fid,prod(nii.hdr.dim(2:nii.hdr.dim(1)+1)),[nii.datatype '=>' nii.datatype]);
fclose(fid);

nii.img=reshape(nii.img,nii.hdr.dim(2:nii.hdr.dim(1)+1));

