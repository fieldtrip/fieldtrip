function spike = ft_datatype_spikeraw(spike, varargin)

% FT_DATATYPE_SPIKERAW describes the FieldTrip MATLAB structure for raw spike data
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
% (2007) Introduced the raw spike data structure.
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
      
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for spike datatype', version);
end



