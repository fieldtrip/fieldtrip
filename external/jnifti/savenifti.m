function bytestream=savenifti(img, filename, varargin)
%
%    savenifti(img, filename)
%        or
%    savenifti(img, filename, rawhdr)
%    savenifti(img, filename, 'nifti2')
%    bytestream=savenifti(img)
%
%    Write an image to a NIfTI (*.nii) or compressed NIfTI file (.nii.gz)
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        img: this is a numerical array to be stored in the NIfTI file
%        filename: output file name, can have a suffix of '.nii' or '.nii.gz'
%                 if a .gz suffix is used, this function needs the JSONLab 
%                 (http://gitlab.com/fangq/jsonlab) and ZMat (http://gitlab.com/fangq/zmat)
%                 to perform the compression.
%        rawhdr (optional): a struct, as a pre-created/loaded NIfTI header data structure
%                if rawhdr is 'nifti1' or 'nifti2', this function calls 
%                nifticreate to create a default header.
%    output:
%        bytestream (optional): the output file byte stream. it only returns this output if
%                no filename is given. 
%
%    example:
%        a=single(rand(10,20,30));
%        savenifti(a,'randnii.nii');
%        savenifti(a,'randnii2.nii.gz','nifti2'); % needs zmat
%
%
%    this file is part of JNIfTI specification: https://github.com/fangq/jnifti
%
%    License: Apache 2.0, see https://github.com/fangq/jnifti for details
%


if(~isempty(varargin))
    if(isstruct(varargin{1}))
        header=varargin{1};
    elseif(ischar(varargin{1}))
        header=nifticreate(img,varargin{1});
    end
else
    header=nifticreate(img);
end

names=fieldnames(header);
buf=[];
for i=1:length(names)
    buf=[buf,typecast(header.(names{i}),'uint8')];
end

if(length(buf)~=352 && length(buf)~=544)
    error('incorrect nifti-1/2 header %d',length(buf));
end

buf=[buf,typecast(img(:)','uint8')];

if(nargout>1 && nargin<2)
    bytestream=buf;
    return;
end

if(regexp(filename,'\.[Gg][Zz]$'))
    buf=gzipencode(buf);
end

fid=fopen(filename,'wb');
if(fid==0)
    error('can not write to the specified file');
end
fwrite(fid,buf);
fclose(fid);
