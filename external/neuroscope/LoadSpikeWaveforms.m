function waveforms = LoadSpikeWaveforms(filename,nChannels,nSamples,list)

%LoadSpikeWaveforms - Read spike waveforms from disk.
%
%  USAGE
%
%    waveforms = LoadSpikeWaveforms(filename,nChannels,nSamples,list)
%
%    filename            spike waveform file name
%    nChannels           number of channels in electrode group
%    nSamples            number of samples per waveform
%    list                optional list of spikes (from 1 to N) to load
%
%  OUTPUT
%
%    waveforms           3D array (spike #,channel,sample) of waveforms
%
%  SEE
%
%    See also GetSpikeWaveforms, PlotSpikeWaveforms.

% Copyright (C) 2004-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help LoadSpikeWaveforms">LoadSpikeWaveforms</a>'' for details).');
end

[path,basename,extension] = fileparts(filename);
if isempty(path), path = '.'; end

electrodeGroupStr = extension(2:end);
electrodeGroup = str2num(electrodeGroupStr);
[~,basename] = fileparts(basename);

% Load .spk file
filename = [path '/' basename '.spk.' electrodeGroupStr];
if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
if nargin < 4,
	waveforms = LoadBinary(filename,'nChannels',nChannels);
else
	waveforms = LoadBinary(filename,'nChannels',nChannels,'nRecords',nSamples*ones(size(list)),'offset',nSamples*(list-1));	
end
waveforms = reshape(waveforms',nChannels,nSamples,[]);
waveforms = permute(waveforms,[3 1 2]); % rearrange: spike #, channel, sample
