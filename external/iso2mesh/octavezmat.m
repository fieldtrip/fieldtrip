function varargout = octavezmat(data, iscompress, zipmethod)
%
% output = octavezmat(input, iscompress, zipmethod)
%    or
% [output, info] = octavezmat(input, iscompress, zipmethod)
% unzipdata = octavezmat(zipdata, info)
%
% Compress or decompress zlib and gzip memory buffers using zip/unzip/gzip/gunzip on Octave
% in case ZMat toolbox (http://github.com/NeuroJSON/zmat) was not installed (ZMat is much faster)
%
% Author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%      input: the input data (can be either compressed or before compression),
%             can be a string, a numerical vector or array
%      iscompress: (optional) if iscompress is 1, zmat compresses/encodes the input,
%             if 0, it decompresses/decodes the input. Default value is 1.
%
%             if iscompress is set to a negative integer, (-iscompress) specifies
%             the compression level. For zlib/gzip, default level is 6 (1-9); for
%             lzma/lzip, default level is 5 (1-9); for lz4hc, default level is 8 (1-16).
%             the default compression level is used if iscompress is set to 1.
%
%             zmat removes the trailing newline when iscompress=2 and method='base64'
%             all newlines are kept when iscompress=3 and method='base64'
%
%             if one defines iscompress as the info struct (2nd output of zmat), zmat
%             will perform a decoding/decompression operation and recover the original
%             input using the info stored in the info structure.
%      method: (optional) compression method, only the below two methods are supported
%             'zlib': zlib/zip based data compression (default)
%             'gzip': gzip formatted data compression
%
% output:
%      output: the decompressed byte stream stored in a uint8 vector; if info is
%            given, output will restore the original data's type and dimensions
%      info: (optional) a struct storing additional info regarding the input data, may have
%            'type': the class of the input array
%            'size': the dimensions of the input array
%            'byte': the number of bytes per element in the input array
%            'method': a copy of the 3rd input indicating the encoding method
%            'status': the zlib/lzma/lz4 compression/decompression function return value,
%                    including potential error codes; see documentation of the respective
%                    libraries for details
%            'level': a copy of the iscompress flag; if non-zero, specifying compression
%                    level, see above
%
% examples:
%      NO_ZMAT=1  % by setting this flag to 1 in the caller or global workspace, octavezmat won't warn zmat is missing
%      [ss,info]=octavezmat(ones(10))
%      orig=octavezmat(ss,info)
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

nowarning = getvarfrom({'caller', 'base'}, 'NO_ZMAT');

if (isempty(nowarning) || nowarning == 0)
    warning('You are recommended to install ZMat (http://github.com/NeuroJSON/zmat) get much faster speed in Octave');
end

if (nargin < 1)
    fprintf(1, 'Format: output = octavezmat(data, iscompress, zipmethod)\n');
    return
end

if (nargin < 2)
    iscompress = 1;
end

if (nargin < 3)
    zipmethod = 'zlib';
end

if (isstruct(iscompress))
    inputinfo = iscompress;
    iscompress = 0;
end

if (~(ischar(data) || islogical(data) || (isnumeric(data) && isreal(data))))
    error('input must be a char, non-complex numeric or logical vector or N-D array');
end

if (ischar(data))
    data = uint8(data);
end

fname = tempname;
tmpfile = fname;
outputfile = fname;

suff = struct('zlib', '.zip', 'gzip', '.gz');

if (~isfield(suff, zipmethod))
    error('zipmethod is not supported');
end

if (~iscompress)
    tmpfile = [fname suff.(zipmethod)];
end

fd = fopen(tmpfile, 'wb');
if (~fd)
    error('unable to create temporary file');
end

fwrite(fd, typecast(data(:), 'uint8'), 'uint8');
fclose(fd);

if (iscompress)
    outputfile = [fname suff.(zipmethod)];
end

if (~iscompress)
    if (strcmp(zipmethod, 'zlib'))
        outputfile = unzip(tmpfile, tempdir);
        if ((exist('OCTAVE_VERSION', 'builtin') ~= 0))
            outputfile = [tempdir filesep outputfile{1}];
        else
            outputfile = outputfile{1};
        end
    elseif (strcmp(zipmethod, 'gzip'))
        gunzip(tmpfile);
    end
else
    if (strcmp(zipmethod, 'zlib'))
        zip(outputfile, tmpfile);
    elseif (strcmp(zipmethod, 'gzip'))
        gzip(tmpfile);
    end
end

if (exist(tmpfile, 'file'))
    delete(tmpfile);
end

fd = fopen(outputfile, 'rb');
if (~fd)
    error('failed to unzip buffer');
end
varargout{1} = fread(fd, [1 inf], 'uint8=>uint8');
fclose(fd);

if (exist(outputfile, 'file'))
    delete(outputfile);
end

if (nargout > 1)
    varargout{2} = struct('type', class(data), 'size', size(data), 'method', zipmethod, 'status', 0, 'level', iscompress);
end

if (exist('inputinfo', 'var') && isfield(inputinfo, 'type'))
    if (strcmp(inputinfo.type, 'logical'))
        varargout{1} = logical(varargout{1});
    else
        varargout{1} = typecast(varargout{1}, inputinfo.type);
    end
    varargout{1} = reshape(varargout{1}, inputinfo.size);
end
