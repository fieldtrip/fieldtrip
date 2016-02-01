function amplitudes = LoadSpikeAmplitudes(filename,nChannels,nSamples,peak,rate)

%LoadSpikeAmplitudes - Read spike times, cluster IDs and amplitudes from disk.
%
%  USAGE
%
%    amplitudes = LoadSpikeAmplitudes(filename,nChannels,nSamples)
%
%    filename            spike waveform file name
%    nChannels           number of channels in electrode group
%    nSamples            number of samples per waveform
%    peak                peak sample (1..N)
%    rate                sampling rate
%
%  SEE
%
%    See also GetSpikeAmplitudes, LoadSpikeWaveforms.

% Copyright (C) 2004-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help LoadSpikeAmplitudes">LoadSpikeAmplitudes</a>'' for details).');
end

[path,basename,extension] = fileparts(filename);
if isempty(path), path = '.'; end

electrodeGroupStr = extension(2:end);
electrodeGroup = str2num(electrodeGroupStr);
[~,basename] = fileparts(basename);

% Times and cluster IDs
times = LoadSpikeTimes(filename,rate);

% Load .spk file
filename = [path '/' basename '.spk.' electrodeGroupStr];
if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
amp = LoadBinary(filename,'offset',peak-1,'skip',nSamples-1,'nChannels',nChannels);

amplitudes = [times amp];
