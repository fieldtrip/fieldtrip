function [timelock] = ft_timelockbaseline(cfg, timelock);

% FT_TIMELOCKBASELINE performs baseline correction for ERF and ERP data
%
% Use as
%    [timelock] = ft_timelockbaseline(cfg, timelock)
% where the timelock data comes from FT_TIMELOCKANALYSIS and the
% configuration should contain
%   cfg.baseline     = [begin end] (default = 'no')
%   cfg.channel      = cell-array, see FT_CHANNELSELECTION
%
% See also FT_TIMELOCKANALYSIS, FT_FREQBASELINE
%
% Undocumented local options:
%   cfg.blcwindow
%   cfg.previous
%   cfg.version
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%   cfg.outputfile = one can specify output as file to save to disk

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

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'baseline'),   cfg.baseline    = 'no';   end
if ~isfield(cfg, 'inputfile'),  cfg.inputfile   = [];     end
if ~isfield(cfg, 'outputfile'), cfg.outputfile  = [];     end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    timelock = loadvar(cfg.inputfile, 'data');
    hasdata = true;
  end
end

% check if the input data is valid for this function
timelock = checkdata(timelock, 'datatype', 'timelock', 'feedback', 'yes');

% the cfg.blc/blcwindow options are used in preprocessing and in
% ft_timelockanalysis (i.e. in private/preproc), hence make sure that
% they can also be used here for consistency
if isfield(cfg, 'baseline') && (isfield(cfg, 'blc') || isfield(cfg, 'blcwindow'))
  error('conflicting configuration options, you should use cfg.baseline');
elseif isfield(cfg, 'blc') && strcmp(cfg.blc, 'no')
  cfg.baseline = 'no';
  cfg = rmfield(cfg, 'blc');
  cfg = rmfield(cfg, 'blcwindow');
elseif isfield(cfg, 'blc') && strcmp(cfg.blc, 'yes')
  cfg.baseline = cfg.blcwindow;
  cfg = rmfield(cfg, 'blc');
  cfg = rmfield(cfg, 'blcwindow');
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

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';

if hasdata && isfield(timelock, 'cfg')
  % remember the configuration details of the input data
  cfg.previous = timelock.cfg;
end
% remember the exact configuration details in the output
timelock.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', timelock); % use the variable name "data" in the output file
end
