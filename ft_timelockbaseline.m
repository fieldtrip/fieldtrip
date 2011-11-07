function [timelock] = ft_timelockbaseline(cfg, timelock)

% FT_TIMELOCKBASELINE performs baseline correction for ERF and ERP data
%
% Use as
%    [timelock] = ft_timelockbaseline(cfg, timelock)
% where the timelock data comes from FT_TIMELOCKANALYSIS and the
% configuration should contain
%   cfg.baseline     = [begin end] (default = 'no')
%   cfg.channel      = cell-array, see FT_CHANNELSELECTION
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_TIMELOCKANALYSIS, FT_FREQBASELINE

% Undocumented local options:
%   cfg.baselinewindow
%   cfg.previous
%   cfg.version

% Copyright (C) 2006, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar timelock

% check if the input data is valid for this function
timelock = ft_checkdata(timelock, 'datatype', 'timelock', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed', {'blcwindow', 'baselinewindow'});

% set the defaults
if ~isfield(cfg, 'baseline'),   cfg.baseline    = 'no';   end

% the cfg.blc/blcwindow options are used in preprocessing and in
% ft_timelockanalysis (i.e. in private/preproc), hence make sure that
% they can also be used here for consistency
if isfield(cfg, 'baseline') && (isfield(cfg, 'demean') || isfield(cfg, 'baselinewindow'))
  error('conflicting configuration options, you should use cfg.baseline');
elseif isfield(cfg, 'demean') && strcmp(cfg.demean, 'no')
  cfg.baseline = 'no';
  cfg = rmfield(cfg, 'demean');
  cfg = rmfield(cfg, 'baselinewindow');
elseif isfield(cfg, 'demean') && strcmp(cfg.demean, 'yes')
  cfg.baseline = cfg.baselinewindow;
  cfg = rmfield(cfg, 'demean');
  cfg = rmfield(cfg, 'baselinewindow');
end

if ischar(cfg.baseline)
  if strcmp(cfg.baseline, 'yes')
    % do correction on the whole time interval
    cfg.baseline = [-inf inf];
  elseif strcmp(cfg.baseline, 'all')
    % do correction on the whole time interval
    cfg.baseline = [-inf inf];
  end
end

if ~(ischar(cfg.baseline) && strcmp(cfg.baseline, 'no'))
  % determine the time interval on which to apply baseline correction
  tbeg = nearest(timelock.time, cfg.baseline(1));
  tend = nearest(timelock.time, cfg.baseline(2));
  % update the configuration
  cfg.baseline(1) = timelock.time(tbeg);
  cfg.baseline(2) = timelock.time(tend);
  
  if isfield(cfg, 'channel')
    % only apply on selected channels
    cfg.channel = ft_channelselection(cfg.channel, timelock.label);
    chansel = match_str(timelock.label, cfg.channel);
    timelock.avg(chansel,:) = ft_preproc_baselinecorrect(timelock.avg(chansel,:), tbeg, tend);
  else
    % apply on all channels
    timelock.avg = ft_preproc_baselinecorrect(timelock.avg, tbeg, tend);
  end
  
  if strcmp(timelock.dimord, 'rpt_chan_time')
    fprintf('applying baseline correction on each individual trial\n');
    ntrial = size(timelock.trial,1);
    if isfield(cfg, 'channel')
      % only apply on selected channels
      for i=1:ntrial
        timelock.trial(i,chansel,:) = ft_preproc_baselinecorrect(shiftdim(timelock.trial(i,chansel,:),1), tbeg, tend);
      end
    else
      % apply on all channels
      for i=1:ntrial
        timelock.trial(i,:,:) = ft_preproc_baselinecorrect(shiftdim(timelock.trial(i,:,:),1), tbeg, tend);
      end
    end
  elseif strcmp(timelock.dimord, 'subj_chan_time')
    fprintf('applying baseline correction on each individual subject\n');
    nsubj = size(timelock.individual,1);
    if isfield(cfg, 'channel')
      % only apply on selected channels
      for i=1:nsubj
        timelock.individual(i,chansel,:) = ft_preproc_baselinecorrect(shiftdim(timelock.individual(i,chansel,:),1), tbeg, tend);
      end
    else
      % apply on all channels
      for i=1:nsubj
        timelock.individual(i,:,:) = ft_preproc_baselinecorrect(shiftdim(timelock.individual(i,:,:),1), tbeg, tend);
      end
    end
  end
  
  if isfield(timelock, 'var')
    fprintf('baseline correction invalidates previous variance estimate, removing var\n');
    timelock = rmfield(timelock, 'var');
  end
  
  if isfield(timelock, 'cov')
    fprintf('baseline correction invalidates previous covariance estimate, removing cov\n');
    timelock = rmfield(timelock, 'cov');
  end
  
end % ~strcmp(cfg.baseline, 'no')

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous timelock
ft_postamble history timelock
ft_postamble savevar timelock

