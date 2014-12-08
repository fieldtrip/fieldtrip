function parameters = LoadParameters(filename)

%LoadParameters - Load parameters from an XML file.
%
%  USAGE
%
%    parameters = LoadParameters(filename)
%
%    filename            parameter file name
%
% This is a slightly modified version from the original: it also contains
% calibration information in the output, added 20141008


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
	parameters.spikeGroups.groups = {};
	parameters.spikeGroups.nGroups = 0;
end

parameters.nBits          = str2double(p.acquisitionSystem.nBits);
parameters.nChannels      = str2double(p.acquisitionSystem.nChannels);
parameters.voltageRange   = str2double(p.acquisitionSystem.voltageRange);
parameters.amplification  = str2double(p.acquisitionSystem.amplification);
parameters.offset         = str2double(p.acquisitionSystem.offset);
parameters.rates.lfp      = str2double(p.fieldPotentials.lfpSamplingRate);
parameters.rates.wideband = str2double(p.acquisitionSystem.samplingRate);
parameters.xmltree        = p;
try
	parameters.rates.video = str2double(p.video.samplingRate);
	parameters.maxX        = str2double(p.video.width);
	parameters.maxY        = str2double(p.video.height);
catch
	parameters.rates.video = 0;
	parameters.maxX        = 0;
	parameters.maxY        = 0;
	disp('... warning: missing video parameters (set to zero)');
end
