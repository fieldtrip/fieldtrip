function outstruct=memmapstream(bytes, format, varargin)
%
%    outstruct=memmapstream(bytes, format)
%
%    Map a byte-array (in char array or uint8/int8 array) into a structure
%    using a dictionary (format is compatible with memmapfile in MATLAB)
%
%    This function is compatible with both MATLAB and GNU Octave. 
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        bytes: a char, int8 or uint8 vector or array
%        format: a 3-column cell array in the format compatible with the
%              'Format' parameter of memmapfile in MATLAB. It has the
%              following structure
%
%             column 1: data type string, it can be one of the following
%                'int8','int16','int32','int64',
%                'uint8','uint16','uint32','uint64',
%                'single','double','logical'
%             column 2: an integer vector denoting the size of the data
%             column 3: a string denoting the fieldname in the output struct
%
%             For example format={'int8',[1,8],'key'; 'float',[1,1],'value'}
%             reads the first 8 bytes from 'bytes' as the first subfield
%             'key' and the following 4 bytes as the floating point 'value'
%             subfield.
%
%    output:
%        outstruct: a structure containing the required field
%
%    example:
%        bytestream=['Andy' 5 'JT'];
%        format={'uint8', [1,4], 'name',
%              'uint8', [1,1], 'age',
%              'uint8', [1,2], 'school'};
%        data=memmapstream(bytestream,format);
%
%    this file is part of JNIfTI specification: https://github.com/fangq/jnifti
%
%    License: Apache 2.0, see https://github.com/fangq/jnifti for details
%

if(nargin<2)
   error('must provide bytes and format as inputs');
end

if(~ischar(bytes) && ~isa(bytes,'int8') && ~isa(bytes,'uint8') || isempty(bytes))
   error('first input, bytes, must be a char-array or uint8/int8 vector');
end

if(~iscell(format) || size(format,2)<3 || size(format,1)==0 || ~ischar(format{1,1}))
   error('second input, format, must be a 3-column cell array, in a format described by the memmapfile Format field.');
end

bytes=bytes(:)';

datatype=struct('int8',1,'int16',2,'int32',4,'int64',8,'uint8',1,'uint16',2,'uint32',4,'uint64',8,'single',4,'double',8);

opt=varargin2struct(varargin{:});
opt.usemap=jsonopt('usemap',0,opt) && exist('containers.Map');

if(opt.usemap)
    outstruct=containers.Map();
else
    outstruct=struct();
end
len=1;
for i=1:size(format,1)
    bytelen=datatype.(format{i,1})*prod(format{i,2});
    if(opt.usemap)
        outstruct(format{i,3})=reshape(typecast(uint8(bytes(len:bytelen+len-1)),format{i,1}),format{i,2});
    else
        outstruct.(format{i,3})=reshape(typecast(uint8(bytes(len:bytelen+len-1)),format{i,1}),format{i,2});
    end
    len=len+bytelen;
    if(len>length(bytes))
        break;
    end
end
