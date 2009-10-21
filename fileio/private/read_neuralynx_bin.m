function [dat] = read_neuralynx_bin(filename, begsample, endsample);

% READ_NEURALYNX_BIN
%
% Use as
%   hdr = read_neuralynx_bin(filename)
% or
%   dat = read_neuralynx_bin(filename, begsample, endsample)
%
% This  is not a formal Neuralynx file format, but at the
% F.C. Donders Centre we use it in conjunction with Neuralynx,
% SPIKESPLITTING and SPIKEDOWNSAMPLE.
%
% The first version of this file format contained in the first 8 bytes the
% channel label as string. Subsequently it contained 32 bit integer values.
%
% The second version of this file format starts with 8 bytes describing (as
% a space-padded string) the data type. The channel label is contained in
% the filename as dataset.chanlabel.bin.
%
% The third version of this file format starts with 7 bytes describing (as
% a zero-padded string) the data type, followed by the 8th byte which
% describes the downscaling for the 8 and 16 bit integer representations.
% The downscaling itself is represented as uint8 and should be interpreted as
% the number of bits to shift. The channel label is contained in the
% filename as dataset.chanlabel.bin.

% Copyright (C) 2007-2008, Robert Oostenveld
%
% $Log: read_neuralynx_bin.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.5  2008/09/30 08:01:04  roboos
% replaced all fread(char=>char) into uint8=>char to ensure that the
% chars are read as 8 bits and not as extended 16 bit characters. The
% 16 bit handling causes problems on some internationalized OS/Matlab
% combinations.
%
% the help of fread specifies "If the precision is 'char' or 'char*1', MATLAB
% reads characters using the encoding scheme associated with the file.
% See FOPEN for more information".
%
% Revision 1.4  2008/07/01 13:00:58  roboos
% optionally read the 16384 byte ascii header instead of hard-coded defaults
%
% Revision 1.3  2008/05/27 13:04:58  roboos
% switched to the 3rd version of the file format, which includes the downscale/calibration value to recover from int32->int16 compression
% added explicit support for version 1, 2 and 3 of the fileformat
% added original details to the output header
%
% Revision 1.2  2007/12/17 16:23:44  roboos
% fixed bug in jumping to correct begin sample
% added support for determining channel name from filename like this "dataset.chanlabel.bin"
% added support for old splitted dma files, which have an 8 byte header with the channel label and are always int32
%
% Revision 1.1  2007/12/12 16:28:42  roboos
% first implementation
%

needhdr = (nargin==1);
needdat = (nargin>=2);

% this is used for backward compatibility
oldformat = false;

% the first 8 bytes contain the header
fid    = fopen(filename, 'rb', 'ieee-le');
magic  = fread(fid, 8, 'uint8=>char')';

% the header describes the format of the subsequent samples
subtype = [];
if     strncmp(magic, 'uint8',   length('uint8'))
  format = 'uint8';
  samplesize = 1;
elseif strncmp(magic, 'int8',    length('int8'))
  format = 'int8';
  samplesize = 1;
elseif strncmp(magic, 'uint16',  length('uint16'))
  format = 'uint16';
  samplesize = 2;
elseif strncmp(magic, 'int16',   length('int16'))
  format = 'int16';
  samplesize = 2;
elseif strncmp(magic, 'uint32',  length('uint32'))
  format = 'uint32';
  samplesize = 4;
elseif strncmp(magic, 'int32',   length('int32'))
  format = 'int32';
  samplesize = 4;
elseif strncmp(magic, 'uint64',  length('uint64'))
  format = 'uint64';
  samplesize = 8;
elseif strncmp(magic, 'int64',   length('int64'))
  format = 'int64';
  samplesize = 8;
elseif strncmp(magic, 'float32', length('float32'))
  format = 'float32';
  samplesize = 4;
elseif strncmp(magic, 'float64', length('float64'))
  format = 'float64';
  samplesize = 8;
else
  warning('could not detect sample format, assuming file format subtype 1 with ''int32''');
  subtype    = 1; % the file format is version 1
  format     = 'int32';
  samplesize = 4;
end

% determine whether the file format is version 2 or 3
if isempty(subtype)
  if all(magic((length(format)+1):end)==' ')
    subtype = 2;
  else
    subtype = 3;
  end
end

% determine the channel name
switch subtype
  case 1
    % the first 8 bytes of the file contain the channel label (padded with spaces)
    label = strtrim(magic);
  case {2, 3}
    % the filename is formatted like "dataset.chanlabel.bin"
    [p, f, x1] = fileparts(filename);
    [p, f, x2] = fileparts(f);
    if isempty(x2)
      warning('could not determine channel label');
      label = 'unknown';
    else
      label = x2(2:end);
    end
    clear p f x1 x2
  otherwise
    error('unknown file format subtype');
end

% determine the downscale factor, i.e. the number of bits that the integer representation has to be shifted back to the left
switch subtype
  case 1
    % these never contained a multiplication factor but always corresponded
    % to the lowest N bits of the original 32 bit integer
    downscale = 0;
  case 2
    % these might contain a multiplication factor but that factor cannot be retrieved from the file
    warning('downscale factor is unknown for ''%s'', assuming that no downscaling was applied', filename);
    downscale = 0;
  case 3
    downscale = double(magic(8));
  otherwise
    error('unknown file format subtype');
end

[p1, f1, x1] = fileparts(filename);
[p2, f2, x2] = fileparts(f1);
headerfile = fullfile(p1, [f2, '.txt']);
if exist(headerfile, 'file')
  orig = neuralynx_getheader(headerfile);
  % construct the header from the accompanying text file
  hdr             = [];
  hdr.Fs          = orig.SamplingFrequency;
  hdr.nChans      = 1;
  hdr.nSamples    = (filesize(filename)-8)/samplesize;
  hdr.nSamplesPre = 0;
  hdr.nTrials     = 1;
  hdr.label       = {label};
else
  % construct the header from the hard-coded defaults
  hdr             = [];
  hdr.Fs          = 32556;
  hdr.nChans      = 1;
  hdr.nSamples    = (filesize(filename)-8)/samplesize;
  hdr.nSamplesPre = 0;
  hdr.nTrials     = 1;
  hdr.label       = {label};
end

if ~needdat
  % also return the file details
  hdr.orig.subtype    = subtype;
  hdr.orig.magic      = magic;
  hdr.orig.format     = format;
  hdr.orig.downscale  = downscale;
  % return only the header details
  dat = hdr;

else
  % read and return the data
  if begsample<1
    begsample = 1;
  end

  if isinf(endsample)
    endsample = hdr.nSamples;
  end

  fseek(fid, 8+(begsample-1)*samplesize, 'bof');   % skip to the beginning of the interesting data
  format = sprintf('%s=>%s', format, format);
  dat = fread(fid, [1 endsample-begsample+1], format);
  if downscale>1
    % the data was downscaled with 2^N, i.e. shifted N bits to the right in case of integer representations
    % now it should be upscaled again with the same amount
    dat = dat.*(2^downscale);  
  end
  if length(dat)<(endsample-begsample+1)
    error('could not read the requested data');
  end
end % needdat

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the file size in bytes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [siz] = filesize(filename)
l = dir(filename);
if l.isdir
  error(sprintf('"%s" is not a file', filename));
end
siz = l.bytes;

