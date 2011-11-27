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
% The spike format is created by using FT_SPIKE_MAKETRIALS which takes a
% SPIKERAW structure as input.
% An example spike data structure is:
%
%         label: {'chan1'  'chan2'  'chan3'}                                  the channel labels
%         time: {[50481x1 double] [50562x1 double] [50537x1 double]}          the time in the trial for each spike (in seconds)
%         trial: {[50481x1 double] [50562x1 double] [50537x1 double]}         the trial in which each spike was observed
%         trialtime: [100x2 double]                                           the begin- and end-time of each trial (in seconds)
%         waveform: {[32x50481 double] [32x50562 double] [32x50537 double]}   spike waveform, here described with 32 samples
%         cfg: [1x1 struct]                                                   the configuration used by the function that generated this data structure
%
% Required fields:
%   - label, trial, time, trialtime
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
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

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
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2007';
end

switch version
  case {'2010','2007'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % there are no changes required
        
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for spike datatype', version);
end



