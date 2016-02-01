function data = LoadBinary(filename,varargin)

%LoadBinary - Load data from a multiplexed binary file.
%
%  Reading a subset of the data can be done in two different manners: either
%  by specifying start time and duration (more intuitive), or by indicating
%  the position and size of the subset in terms of number of records (more
%  accurate) - a 'record' is a chunk of data containing one sample for each
%  channel.
%
%  LoadBinary can also deal with lists of start times and durations (or
%  offsets and number of records).
%
%  USAGE
%
%    data = LoadBinary(filename,<options>)
%
%    filename       file to read
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz, default = 20kHz)
%     'start'       position to start reading (in s, default = 0)
%     'duration'    duration to read (in s, default = Inf)
%     'offset'      position to start reading (in records, default = 0)
%     'nRecords'    number of records to read (default = Inf)
%     'samples'     same as above (for backward compatibility reasons)
%     'nChannels'   number of data channels in the file (default = 1)
%     'channels'    channels to read (default = all)
%     'precision'   sample precision (default = 'int16')
%     'skip'        number of records to skip after each record is read
%                   (default = 0)
%    =========================================================================

% Copyright (C) 2004-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
nChannels = 1;
precision = 'int16';
skip = 0;
frequency = 20000;
channels = [];
start = 0;
duration = Inf;
offset = 0;
nRecords = Inf;
time = false;
records = false;

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+3) ' is not a property (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'start',
			start = varargin{i+1};
			if ~isdvector(start),
				error('Incorrect value for property ''start'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			if start < 0, start = 0; end
			time = true;
		case 'duration',
			duration = varargin{i+1};
			if ~isdvector(duration,'>=0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			time = true;
		case 'offset',
			offset = varargin{i+1};
			if ~isivector(offset),
				error('Incorrect value for property ''offset'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			if offset < 0, offset = 0; end
			records = true;
		case {'nrecords','samples'},
			nRecords = varargin{i+1};
			if ~isivector(nRecords,'>=0'),
				error('Incorrect value for property ''nRecords'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			if length(nRecords) > 1 && any(isinf(nRecords(1:end-1))),
				error('Incorrect value for property ''nRecords'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			records = true;
		case 'nchannels',
			nChannels = varargin{i+1};
			if ~isiscalar(nChannels,'>0'),
				error('Incorrect value for property ''nChannels'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'channels',
			channels = varargin{i+1};
			if ~isivector(channels,'>=0'),
				error('Incorrect value for property ''channels'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'precision',
			precision = varargin{i+1};
			if ~isstring(precision),
				error('Incorrect value for property ''precision'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'skip',
			skip = varargin{i+1};
			if ~isiscalar(skip,'>=0'),
				error('Incorrect value for property ''skip'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).']);
	end
end

% Either start+duration, or offset+size
if time && records,
	error(['Data subset can be specified either in time or in records, but not both (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).']);
end

% By default, load all channels
if isempty(channels),
	channels = 1:nChannels;
end

% Check consistency between channel IDs and number of channels
if any(channels>nChannels),
	error('Cannot load specified channels (listed channel IDs inconsistent with total number of channels).');
end

% Open file
if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
f = fopen(filename,'r');
if f == -1,
	error(['Cannot read ' filename ' (insufficient access rights?).']);
end

% Size of one data point (in bytes)
sampleSize = 0;
switch precision,
	case {'uchar','unsigned char','schar','signed char','int8','integer*1','uint8','integer*1'},
		sampleSize = 1;
	case {'int16','integer*2','uint16','integer*2'},
		sampleSize = 2;
	case {'int32','integer*4','uint32','integer*4','single','real*4','float32','real*4'},
		sampleSize = 4;
	case {'int64','integer*8','uint64','integer*8','double','real*8','float64','real*8'},
		sampleSize = 8;
end

% Position and number of records of the data subset
if time,
	if length(duration) == 1,
		duration = repmat(duration,size(start,1),1);
	elseif length(duration) ~= length(start),
		error('Start and duration lists have different lengths (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
	end
	dataOffset = floor(start*frequency)*nChannels*sampleSize;
	nRecords = floor(duration*frequency);
else
	if length(nRecords) == 1,
		nRecords = repmat(nRecords,size(offset,1),1);
	elseif length(nRecords) ~= length(offset),
		error('Offset and number of records lists have different lengths (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
	end
	dataOffset = offset*nChannels*sampleSize;
end

% Determine total number of records in file
fileStart = ftell(f);
status = fseek(f,0,'eof');
if status ~= 0,
	fclose(f);
	error('Error reading the data file (possible reasons include trying to read past the end of the file).');
end
fileStop = ftell(f);

% Last number of records may be infinite, compute explicit value
if isinf(nRecords(end)),
	status = fseek(f,dataOffset(end),'bof');
	if status ~= 0,
		fclose(f);
		error('Error reading the data file (possible reasons include trying to read past the end of the file).');
	end
	lastOffset = ftell(f);
	lastNRecords = floor((fileStop-lastOffset)/nChannels/sampleSize);
	nRecords(end) = lastNRecords;
end

% Preallocate memory
data = zeros(sum(nRecords)/(skip+1),length(channels));

% Loop through list of start+duration or offset+nRecords
i = 1;
for k = 1:length(dataOffset),

	% Position file index for reading
	status = fseek(f,dataOffset(k),'bof');
	fileOffset = ftell(f);
	if status ~= 0,
		fclose(f);
		error('Could not start reading (possible reasons include trying to read past the end of the file).');
	end

	% (floor in case all channels do not have the same number of samples)
	maxNRecords = floor((fileStop-fileOffset)/nChannels/sampleSize);
	if nRecords(k) > maxNRecords, nRecords(k) = maxNRecords; end

	% For large amounts of data, read chunk by chunk
	maxSamplesPerChunk = 10000;
	nSamples = nRecords(k)*nChannels;
	if nSamples <= maxSamplesPerChunk,
		d = LoadChunk(f,nChannels,channels,nRecords(k),precision,skip*sampleSize);
		[m,n] = size(d);
		if m == 0, break; end
		data(i:i+m-1,:) = d;
		i = i+m;
	else
		% Determine chunk duration and number of chunks
		nSamplesPerChunk = floor(maxSamplesPerChunk/nChannels)*nChannels;
		nChunks = floor(nSamples/nSamplesPerChunk)/(skip+1);
		% Read all chunks
		for j = 1:nChunks,
			d = LoadChunk(f,nChannels,channels,nSamplesPerChunk/nChannels,precision,skip*sampleSize);
			[m,n] = size(d);
			if m == 0, break; end
			data(i:i+m-1,:) = d;
			i = i+m;
		end
		% If the data size is not a multiple of the chunk size, read the remainder
		remainder = nSamples - nChunks*nSamplesPerChunk;
		if remainder ~= 0,
			d = LoadChunk(f,nChannels,channels,remainder/nChannels,precision,skip*sampleSize);
			[m,n] = size(d);
			if m == 0, break; end
			data(i:i+m-1,:) = d;
			i = i+m;
		end
	end
end

fclose(f);

% ---------------------------------------------------------------------------------------------------------

function data = LoadChunk(fid,nChannels,channels,nSamples,precision,skip)

if skip ~= 0,
	data = fread(fid,[nChannels nSamples],[int2str(nChannels) '*' precision],skip*nChannels);
else
	data = fread(fid,[nChannels nSamples],precision);
end
data = data';

if isempty(data),
	warning('No data read (trying to read past file end?)');
elseif ~isempty(channels),
	data = data(:,channels);
end
