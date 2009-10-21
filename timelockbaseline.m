function [timelock] = timelockbaseline(cfg, timelock);

% TIMELOCKBASELINE performs baseline correction for ERF and ERP data
%
% Use as
%    [timelock] = timelockbaseline(cfg, timelock)
% where the timelock data comes from TIMELOCKANALYSIS and the
% configuration should contain
%   cfg.baseline     = [begin end] (default = 'no')
%   cfg.channel      = cell-array, see CHANNELSELECTION
%
% See also TIMELOCKANALYSIS, FREQBASELINE

% Undocumented local options:
% cfg.blcwindow
% cfg.previous
% cfg.version

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: timelockbaseline.m,v $
% Revision 1.13  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.12  2008/11/21 10:39:10  sashae
% added call to checkconfig
%
% Revision 1.11  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.10  2008/06/10 16:45:05  sashae
% replaced call to blc function with preproc_baselinecorrect
%
% Revision 1.9  2007/12/18 14:17:29  ingnie
% Changed timelock.subj to timelock.individual (Thanks to Ian)
%
% Revision 1.8  2007/06/15 10:00:42  jansch
% fixed typo thanks to Ian
%
% Revision 1.7  2007/05/01 09:19:04  roboos
% fixed bug in removal of covariance (thanks to Marie)
%
% Revision 1.6  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.5  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.4  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.3  2006/05/09 12:22:41  ingnie
% updated help
%
% Revision 1.2  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.1  2006/03/09 10:34:17  roboos
% ne wimplementation, used for plotting functions (c.f. freqbaseline)
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
timelock = checkdata(timelock, 'datatype', 'timelock', 'feedback', 'yes');

% the cfg.blc/blcwindow options are used in preprocessing and in
% timelockanalysis (i.e. in private/preproc), hence make sure that
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

% set the defaults
if ~isfield(cfg, 'baseline'), cfg.baseline = 'no'; end

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
    cfg.channel = channelselection(cfg.channel, timelock.label);
    chansel = match_str(timelock.label, cfg.channel);
    timelock.avg(chansel,:) = preproc_baselinecorrect(timelock.avg(chansel,:), tbeg, tend);
  else
    % apply on all channels
    timelock.avg = preproc_baselinecorrect(timelock.avg, tbeg, tend);
  end

  if strcmp(timelock.dimord, 'rpt_chan_time')
    fprintf('applying baseline correction on each individual trial\n');
    ntrial = size(timelock.trial,1);
    if isfield(cfg, 'channel')
      % only apply on selected channels
      for i=1:ntrial
        timelock.trial(i,chansel,:) = preproc_baselinecorrect(shiftdim(timelock.trial(i,chansel,:),1), tbeg, tend);
      end
    else
      % apply on all channels
      for i=1:ntrial
        timelock.trial(i,:,:) = preproc_baselinecorrect(shiftdim(timelock.trial(i,:,:),1), tbeg, tend);
      end
    end
  elseif strcmp(timelock.dimord, 'subj_chan_time')
    fprintf('applying baseline correction on each individual subject\n');
    nsubj = size(timelock.individual,1);
    if isfield(cfg, 'channel')
      % only apply on selected channels
      for i=1:nsubj
        timelock.individual(i,chansel,:) = preproc_baselinecorrect(shiftdim(timelock.individual(i,chansel,:),1), tbeg, tend);
      end
    else
      % apply on all channels
      for i=1:nsubj
        timelock.individual(i,:,:) = preproc_baselinecorrect(shiftdim(timelock.individual(i,:,:),1), tbeg, tend);
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
cfg.version.id = '$Id: timelockbaseline.m,v 1.13 2009/01/20 13:01:31 sashae Exp $';
% remember the configuration details of the input data
try, cfg.previous = timelock.cfg; end
% remember the exact configuration details in the output 
timelock.cfg = cfg;

