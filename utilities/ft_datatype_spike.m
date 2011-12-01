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
% The first characteristic of the SPIKE structure is a .label cell array
% with the label of the (single or multi) units.
%
%         label: {'unit1'  'unit2'  'unit3'}
%         
% The main fields of the SPIKE structure then contain specific information
% per spike, for the different fields.
% 
% The SPIKE structure must then either contain a .timestamp (cell array)
% field, having the raw timestamps as recorded by the hardware system, or a
% .time (cell array) field, containing the spike times relative to an experimental
% trigger.
%
% For example, the original content of the .timestamp field can be
%
%         timestamp: {[1x993 uint64]  [1x423 uint64]  [1x3424 uint64]}
%
%
% An additional, optional field that is typically obtained from the raw recording is
% the .waveform field
% This again is a cell-array with the waveforms for every unit and label.
% 
% For example, the content of this field may be
%
%         waveform: {[32x993 double] [32x423 double] [32x3224 double]}
%
% Using FT_SPIKE_REDEFINETRIAL we can make the spike times relative to a
% trigger.
%
%
% The .time field, a cell-array, then contains less spikes, and spike times
% are now relative to the trigger.
%
%         time: {[1x504 double] [1x50 double] [1x101 double]}                   the time in the trial for each spike relative to the trigger (in seconds)
%
% In addition, for every trial we register in which trial the spike was
% recorded, this is specified by the .trial field.
%
%         trial: {[1x504 double] [1x50 double] [1x101 double]}                   the trial in which each spike was observed
%
% To fully reconstruct the structure of the spike-train, it is required
% that the starts and the ends of the trials are taken into account.
%
% This is specified by the .trialtime field, a nTrials x 2 matrix.
%
%         trialtime: [100x2 double]                                           the begin- and end-time of each trial (in seconds)
%
%
% From these trials the spike-train can be fully reconstructed.
%
% For example after running FT_SPIKE_REDEFINETRIALS the obtained structure
% may be:
%         label: {'unit1'  'unit2'  'unit3'}
%         timestamp: {[1x504 double] [1x50 double] [1x101 double]} 
%         time: {[1x504 double] [1x50 double] [1x101 double]}  
%         trial: {[1x504 double] [1x50 double] [1x101 double]}  
%         trialtime: [100x2 double]                
%         waveform: {[32x504 double] [32x50 double] [32x101 double]}
%         cfg
%
% For analysing the spike-LFP phase-locking, the output structure is also a
% SPIKE structure with additional fields .fourierspctrm, .lfplabel, .freq and
% .fourierspctrmdimord.
%
% For example, from the structure above we may obtain
%
%         label: {'unit1'  'unit2'  'unit3'}
%         timestamp: {[1x504 double] [1x50 double] [1x101 double]} 
%         time: {[1x504 double] [1x50 double] [1x101 double]}  
%         trial: {[1x504 double] [1x50 double] [1x101 double]}  
%         trialtime: [100x2 double]                
%         waveform: {[32x504 double] [32x50 double] [32x101 double]}
%         fourierspctrmdimord: {chan}_spike_lfplabel_freq
%         fourierspctrm: {504x2x20, 50x2x20, 101x2x20}
%         freq: [1x20 double]
%         lfplabel: {'lfpchan1', 'lfpchan2'}
%         cfg
%
% Required fields:
%   - label, timestamp or the combination (time,trial,trialtime)
%
% Optional fields:
%   - waveform
%   - fourierspctrm, fourierspctrmdimord,freq,lfplabel      extra outputs from FT_SPIKETRIGGEREDSPECTRUM and FT_SPIKE_TRIGGEREDSPECTRUM
%   - hdr
%   - cfg
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

% Copyright (C) 2011, Robert Oostenveld & Martin Vinck
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
          
    if isfield(spike,'origtrial') && isfield(spike,'origtime')
      warning('The spike datatype format you are using is depreciated. Converting to newer spike format');
      spike.trial = {spike.origtrial};
      spike.time  = {spike.origtime};
      if ~isa(spike.fourierspctrm, 'cell')
        spike.fourierspctrm = {spike.fourierspctrm};
      end
      if ~isfield(spike, 'trialtime')
        % determine from the data itself
        warning('Reconstructing the field trialtime from spike.origtime and spike.origtrial');      
        tmax  = nanmin(spike.trial{1});
        tsmin = nanmin(spike.time{1});
        tsmax = nanmax(spike.time{1});
        spike.trialtime = [tsmin*ones(tmax,1) tsmax*ones(tmax,1)];
      end
      spike.lfplabel = spike.label; % in the old format, these were the lfp channels
      spike.label    = {'unit1'};
      spike.fourierspctrmdimord = '{chan}_spike_lfplabel_freq';      
    end
    
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for spike datatype', version);
end



