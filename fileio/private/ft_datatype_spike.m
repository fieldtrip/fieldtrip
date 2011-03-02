function spike = ft_datatype_spike(spike, varargin)

% FT_DATATYPE_SPIKE describes the FieldTrip MATLAB structure for spike data
%
% Spike data is characterised as a sparse point-process, i.e. each neuronal
% firing is only represented as the time at which the firing happened.
% Optionally, the spike waveform can also be represented. Using the spike
% waveform, the neuronal firing events can be sorted into their single units.
%
% Spike data is obtained using FT_READ_SPIKE to read it from a Plexon,
% Neuralynx or other animal electrophysiology system file format containing
% spikes.
%
% An example spike data structure is:
%
%         label: {'chan1', 'chan2', 'chan3'}                         the channel labels
%     timestamp: {[1x993 uint64]  [1x423 uint64]  [1x3424 uint64]}   timestamp in arbitrary units, depends on the acquisition system
%      waveform: {[32x993 double] [32x433 double] [32x3424 double]}  spike waveform, here described with 32 samples
%           hdr: [1x1 struct]                                        the full header information of the original dataset on disk
%
% The example above contains three spike channels, each with varying numbers of
% detected spikes. The timestamps of the spikes are represented, with their
% waveforms (32 samples per waveform). This type of representation can be seen
% as raw spike data, because there is no reference to the experimental trials.
%
% The spike data representation can also represent spike timepoints in relation
% to the experimental trials, as in this example:
%
%         label: {'chan1'  'chan2'  'chan3'}                           the channel labels
%          time: {[50481x1 double] [50562x1 double] [50537x1 double]}  the time in the trial for each spike (in seconds)
%         trial: {[50481x1 double] [50562x1 double] [50537x1 double]}  the trial in which each spike was observed
%     trialtime: [100x2 double]                                        the begin- and end-time of each trial (in seconds)
%           cfg: [1x1 struct]                                          the configuration used by the function that generated this data structure
%
% Required fields:
%   - label, timestamp
%
% Optional fields:
%   - waveform, hdr, cfg
%
% Deprecated fields:
%   - <unknown>
%
% Obsoleted fields:
%   - <unknown>
%
% Revision history:
%
% (2010/latest) Introduced the time and the trialtime fields.
%
% (2007) Introduced the spike data structure.
%
% See also FT_DATATYPE and FT_DATATYPE_xxx

% Copyright (C) 2011, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$


% get the optional input arguments, which should be specified as key-value pairs
version = keyval('version', varargin); if isempty(version), version = 'latest'; end

if strcmp(version, 'latest')
  version = '2007';
end

switch version
  case '2010'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % there are no changes required
    
  case '2007'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the 2007 format did not contain these fields
    if isfield(spike, 'time')
      spike = rmfield(spike, 'time');
    end
    if isfield(spike, 'trialtime')
      spike = rmfield(spike, 'trialtime');
    end
    
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for spike datatype', version);
end



