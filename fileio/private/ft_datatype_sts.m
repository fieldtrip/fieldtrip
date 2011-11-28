function sts = ft_datatype_spike(sts, varargin)

% FT_DATATYPE_STS describes the FieldTrip MATLAB structure for
% spiketriggered spectrum data
%
% Spike triggered spectrum data is characterised as a sparse point-process, i.e. each neuronal
% firing is only represented as the time at which the firing happened.
% For every spike there is a phase value for every (chan x freq)
% combination
%
% sts data is obtained from FT_SPIKE_TRIGGEREDSPECTRUM or
% FT_SPIKETRIGGEREDSPECTRUM
%
% An example sts data structure is:
%
%         label: {'lfpchan1'  'lfpchan2'}                                     the channel labels
%         time: {[50481x1 double] [50562x1 double] [50537x1 double]}          the time in the trial for each spike (in seconds)
%         trial: {[50481x1 double] [50562x1 double] [50537x1 double]}         the trial in which each spike was observed
%         trialtime: [100x2 double]                                           the begin- and end-time of each trial (in seconds)
%         spikechannel: {'unit1', 'unit2', 'unit3'};                          the labels of the single units        
%         fourierspctrm: {[50481x2x20],[50562x2x20],[50537x2x20]};            the fourier spectra for every spike          
%         freq: [1x20 double];                                                the vector of frequencies
%         cfg: [1x1 struct]                                                   the configuration used by the function that generated this data structure
%
% Required fields:
%   - label, trial, time, trialtime, fourierspctrm, spikechannel, freq
%
% Optional fields:
%   - hdr, cfg
%
% Deprecated fields:
%   - origtrial, origtime
%
% Revision history:
%
% (2011/latest) Made the sts structure similar to the spike structure
%
% (2008): introduced the sts format
%
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2011, Martin Vinck
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
% $Id: ft_datatype_sts.m 4866 2011-11-28 09:27:06Z marvin $


% get the optional input arguments, which should be specified as key-value pairs
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2007';
end

switch version
  case {'2010','2007'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % there are no changes required
   
    if isfield(sts,'origtrial') && isfield(sts,'origtime')
      warning('The sts datatype format you are using is depreciated. Converting to newer sts format');
      sts.trial = {sts.origtrial};
      sts.time  = {sts.origtime};
      if ~isa(sts.fourierspctrm, 'cell')
        sts.fourierspctrm = {sts.fourierspctrm};
      end
      if ~isfield(sts, 'spikechannel')
        sts.spikechannel = {'unit1'};
      end
      if ~isfield(sts, 'trialtime')
        % determine from the data itself
        tmax  = nanmin(sts.trial{1});
        tsmin = nanmin(sts.time{1});
        tsmax = nanmax(sts.time{1});
        sts.trialtime = [tsmin*ones(tmax,1) tsmax*ones(tmax,1)];
      end
    end
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for spike datatype', version);
end



