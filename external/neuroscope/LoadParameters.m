function parameters = LoadParameters(filename)

%LoadParameters - Load parameters from an XML file.
%
%  USAGE
%
%    parameters = LoadParameters(filename)
%
%    filename            parameter file name

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
[pathname,basename] = fileparts(filename);
parameters.session.path = pathname;
parameters.session.name = basename;

if isempty(which('xmltree')),
	error('This function requires the <a href="http://www.artefact.tk/software/matlab/xml/">xmltree</a> toolbox by G. Flandin.');
end

t = xmltree(filename);
p = convert(t);
if ~isempty(p.spikeDetection),
	parameters.spikeGroups.nGroups = length(p.spikeDetection.channelGroups.group);
	if parameters.spikeGroups.nGroups == 1,
		parameters.spikeGroups.nSamples = str2num(p.spikeDetection.channelGroups.group.nSamples);
		parameters.spikeGroups.peakSamples = str2num(p.spikeDetection.channelGroups.group.peakSampleIndex);
		channels = p.spikeDetection.channelGroups.group.channels.channel;
		if isa(channels,'cell'),
			for channel = 1:length(channels),
				parameters.spikeGroups.groups{1}(channel) = str2num(channels{channel});
			end
		else
			parameters.spikeGroups.groups{1} = str2num(channels);
		end
	else
		for group = 1:parameters.spikeGroups.nGroups,
			parameters.spikeGroups.nSamples(group) = str2num(p.spikeDetection.channelGroups.group{group}.nSamples);
			parameters.spikeGroups.peakSamples(group) = str2num(p.spikeDetection.channelGroups.group{group}.peakSampleIndex);
			channels = p.spikeDetection.channelGroups.group{group}.channels.channel;
			if isa(channels,'cell'),
				for channel = 1:length(channels),
					parameters.spikeGroups.groups{group}(channel) = str2num(channels{channel});
				end
			else
				parameters.spikeGroups.groups{group} = str2num(channels);
			end
		end
	end
else
	parameters.spikeGroups.nSamples = [];
	parameters.spikeGroups.peakSampleIndex = [];
	parameters.spikeGroups.groups = {};
	parameters.spikeGroups.nGroups = 0;
end

parameters.nChannels = str2num(p.acquisitionSystem.nChannels);
parameters.nBits = str2num(p.acquisitionSystem.nBits);
parameters.rates.lfp = str2num(p.fieldPotentials.lfpSamplingRate);
parameters.rates.wideband = str2num(p.acquisitionSystem.samplingRate);
try
	parameters.rates.video = str2num(p.video.samplingRate);
	parameters.maxX = str2num(p.video.width);
	parameters.maxY = str2num(p.video.height);
catch
	parameters.rates.video = 0;
	parameters.maxX = 0;
	parameters.maxY = 0;
	disp('... warning: missing video parameters (set to zero)');
end
