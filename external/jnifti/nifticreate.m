function header=nifticreate(img,format)
%
%    header=nifticreate(img)
%       or
%    header=nifticreate(img,format)
%
%    Create a default NIfTI header
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        img: the image data matching the header
%        format: can only be 'nifti1'. can be ignored
%
%    output:
%        header: a struct that is byte-wise compatible with NIfTI-1
%
%    this file is part of JNIfTI specification: https://github.com/fangq/jnifti
%
%    License: Apache 2.0, see https://github.com/fangq/jnifti for details
%

if(nargin<2)
    format='nifti1';
end

datatype=struct('int8',256,'int16',4,'int32',8,'int64',1024,'uint8',2,'uint16',512,'uint32',768,'uint64',1280,'single',16,'double',64);

if(strcmp(format,'nifti1'))
    headerlen=348;
else
    headerlen=540;
end

header=memmapstream(uint8(zeros(1,headerlen+4)),niiformat(format));
header.sizeof_hdr=cast(headerlen, class(header.sizeof_hdr));
header.datatype=cast(datatype.(class(img)), class(header.datatype));
header.dim(1:end)=cast(1,class(header.dim));
header.dim(1:ndims(img)+1)=cast([ndims(img),size(img)], class(header.dim));
header.pixdim(1:ndims(img)+1)=cast(1, class(header.pixdim));
header.vox_offset=cast(headerlen+4, class(header.vox_offset));
if(header.sizeof_hdr==540)
    header.magic(1:3)=cast('ni2',class(header.magic));
else
    header.magic(1:3)=cast('ni1',class(header.magic));
end
header.srow_x(1)=cast(1, class(header.srow_x));
header.srow_y(2)=cast(1, class(header.srow_y));
header.srow_z(3)=cast(1, class(header.srow_z));
header.sform_code=cast(1, class(header.sform_code));
