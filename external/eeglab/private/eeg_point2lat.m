% eeg_point2lat() - convert latency in data points to latency in ms relative
%                   to the time locking. Used in eeglab().
% Usage:
%       >> [newlat ] = eeg_point2lat( lat_array, [], srate);
%       >> [newlat ] = eeg_point2lat( lat_array, epoch_array,...
%                                 srate, timelimits, timeunit);
% Inputs:
%   lat_array   - latency array in data points assuming concatenated
%                 data epochs (see eeglab() event structure)
%   epoch_array - epoch number corresponding to each latency value
%   srate       - data sampling rate in Hz
%   timelimits  - [min max] timelimits in 'timeunit' units (see below)
%   timeunit    - time unit in second. Default is 1 = seconds.
%
% Outputs:
%   newlat      - converted latency values (in 'timeunit' units) for each epoch
%
% Example:
%   tmpevent = EEG.event;
%   eeg_point2lat( [ tmpevent.latency ], [], EEG.srate, [EEG.xmin EEG.xmax]);
%   % returns the latency of all events in second for a continuous
%   % dataset EEG
%
%   eeg_point2lat( [ tmpevent.latency ], [ tmpevent.epoch ], 
%                 EEG.srate, [EEG.xmin EEG.xmax]*1000, 1E-3);
%   % returns the latency of all events in millisecond for a dataset
%   % containing data epochs.
%
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2 Mai 2002
%
% See also: eeg_lat2point(), eeglab(), pop_editieventvals(), pop_loaddat()

% Copyright (C) 2 Mai 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function newlat = eeg_point2lat( lat_array, epoch_array, srate, timewin, timeunit);

if nargin <3
    help eeg_point2lat;
    return;
end
if isempty( epoch_array )
    epoch_array = ones( size(lat_array) );
end
if nargin <4
    timewin = 0;
end
if nargin <5
	timeunit = 1;
end

if length(lat_array) ~= length(epoch_array)
	if length(epoch_array)~= 1
		disp('eeg_point2lat: latency and epoch arrays must have the same length'); return;
	else
		epoch_array = ones(1,length(lat_array))*epoch_array;
	end
end
if length(timewin) ~= 2
    disp('eeg_point2lat: timelimits array must have length 2'); return;
end
if iscell(epoch_array)
	epoch_array = [ epoch_array{:} ];
end
if iscell(lat_array)
	lat_array = [ lat_array{:} ];
end

timewin = timewin*timeunit;

if length(timewin) == 2
    pnts = (timewin(2)-timewin(1))*srate+1;
else
    pnts = 0;
end
newlat  = ((lat_array - (epoch_array-1)*pnts-1)/srate+timewin(1))/timeunit;
newlat = round(newlat*1E9)*1E-9;
