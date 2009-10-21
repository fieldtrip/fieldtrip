function [grandavg] = timelockgrandaverage(cfg, varargin)

% TIMELOCKGRANDAVERAGE computes ERF/ERP average and variance
% over multiple subjects
%
% Use as
%   [grandavg] = timelockgrandaverage(cfg, avg1, avg2, avg3, ...)
%
% where
%   avg1..N are the ERF/ERP averages as obtained from TIMELOCKANALYSIS
%
% and cfg is a configuration structure with
%  cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                       see CHANNELSELECTION for details
%  cfg.latency        = [begin end] in seconds or 'all' (default = 'all')
%  cfg.keepindividual = 'yes' or 'no' (default = 'no')
%  cfg.normalizevar   = 'N' or 'N-1' (default = 'N-1')
%
% See also TIMELOCKANALYSIS, TIMELOCKSTATISTICS

% Copyright (C) 2003-2006, Jens Schwarzbach
%
% $Log: timelockgrandaverage.m,v $
% Revision 1.21  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.20  2008/11/21 12:48:17  sashae
% added call to checkconfig at start and end of function
%
% Revision 1.19  2008/10/21 09:40:05  roboos
% reset min/max latency for variable trial length over subjects
% shift start and end index in case of flipped time axis, potentially introduced for response time locked data
% thanks to Christian & Ines
%
% Revision 1.18  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.17  2007/12/12 07:54:00  roboos
% give warning if discarding grad or elec
%
% Revision 1.16  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.15  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.14  2006/10/04 07:10:08  roboos
% updated documentation
%
% Revision 1.13  2006/08/29 15:05:26  roboos
% fixed a bug in channel selection, the selection was not the intersection of all inputs but the combination of them (which did not match)
%
% Revision 1.12  2006/06/20 16:25:58  ingnie
% updated documentation
%
% Revision 1.11  2006/06/13 14:48:11  ingnie
% updated documentation
%
% Revision 1.10  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.9  2005/06/02 12:12:15  roboos
% changed indentation and whitespace
% changed the dimord *_tim into *_time
% removed nsubjects from output
%
% Revision 1.8  2005/05/17 17:50:39  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.7  2005/04/25 09:43:07  roboos
% fixed bug for different ordering of channels in input data
% replaced squeeze() with reshape()
%
% Revision 1.6  2004/10/13 14:11:30  roboos
% changed cfg.previous, now consistent over functions with try-statement
% and will also work if the input already had a cfg.previous
%
% Revision 1.5  2004/09/16 15:36:51  roboos
% fixed bug in latency selection
%
% Revision 1.4  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details
%
% Revision 1.3  2004/06/30 15:05:53  roboos
% converted ascii file from DOS to UNIX format, no functional code changes
%
% Revision 1.2  2004/06/29 15:42:39  roboos
% added CVS Log feature
% fixed bug in latency selection for cfg.latency=all
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'no');
end

% set the defaults
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';		   end
if ~isfield(cfg, 'keepindividual'), cfg.keepindividual = 'no'; end
if ~isfield(cfg, 'latency'),        cfg.latency = 'all';	     end
if ~isfield(cfg, 'normalizevar'),   cfg.normalizevar = 'N-1';  end

% varargin{1} ... varargin{end} contain the individual ERFs
Nsubj = length(varargin);

if isfield(varargin{1}, 'grad')
  warning('discarding gradiometer information because it cannot be averaged');
end
if isfield(varargin{1}, 'elec')
  warning('discarding electrode information because it cannot be averaged');
end

% replace string latency selection by a timerange based on the range of all subjects
if ischar(cfg.latency) && strcmp(cfg.latency, 'all')
  cfg.latency = [];
  cfg.latency(1) = min(varargin{1}.time);
  cfg.latency(2) = max(varargin{1}.time);
  for s=2:Nsubj
    % reset min/max latency (for variable trial length over subjects)
    if min(varargin{s}.time) > cfg.latency(1)
      cfg.latency(1) = min(varargin{s}.time);
    end
    if max(varargin{s}.time) < cfg.latency(2)
      cfg.latency(2) = max(varargin{s}.time);
    end
  end
end

%SELECT TIME WINDOW
idxs = nearest(varargin{1}.time, min(cfg.latency));
idxe = nearest(varargin{1}.time, max(cfg.latency));
% shift start and end index in case of flipped time axis (potentially introduced for response time locked data)
if idxe < idxs
  ResultsTimeSelectCases = idxe:idxs;
else
  ResultsTimeSelectCases = idxs:idxe;
end
ResultsNTimePoints = length(ResultsTimeSelectCases);
ResultsTime = varargin{1}.time(ResultsTimeSelectCases);

%UPDATE CFG STRUCTURE WITH TIME THAT WAS FINALLY USED
cfg.latency = [ResultsTime(1), ResultsTime(end)];

%DETERMINE WHICH CHANNELS ARE AVAILABLE FOR ALL SUBJECTS
for i=1:Nsubj
  cfg.channel = channelselection(cfg.channel, varargin{i}.label);
end
ResultNChannels = size(cfg.channel, 1);

%REDUCE DATASET TO INTERSECTION OF DESIRED AND AVAILABLE CHANNELS
for i=1:Nsubj
  % select channel indices in this average, sorted according to configuration
  [dum, chansel]    = match_str(cfg.channel, varargin{i}.label);
  varargin{i}.avg   = varargin{i}.avg(chansel,:);
  varargin{i}.label = varargin{i}.label(chansel);
  try, varargin{i}.trial = varargin{i}.trial(chansel,:,:); end
end

%PREALLOCATE
avgmat = zeros(Nsubj, ResultNChannels, ResultsNTimePoints);
%FILL MATRIX, MAY BE DONE MORE EFFECTIVELY WITH DEAL COMMAND
for s = 1:Nsubj
  avgmat(s, :, :) = varargin{s}.avg(:, ResultsTimeSelectCases);
end

%AVERAGE ACROSS SUBJECT DIMENSION
ResultGrandavg = mean(avgmat, 1);
ResultGrandavg = reshape(ResultGrandavg, [ResultNChannels, ResultsNTimePoints]);

%COMPUTE VARIANCE ACROSS SUBJECT DIMENSION
%THIS LOOKS AWKWARD (std.^2) BUT IS FAST DUE TO BUILT IN FUNCTIONS
switch cfg.normalizevar
  case 'N-1'
    sdflag = 0;
  case 'N'
    sdflag = 1;
end
ResultVar = std(avgmat, sdflag, 1).^2;
ResultVar = reshape(ResultVar, [ResultNChannels, ResultsNTimePoints]);

%--------------------------------------------
% % collect the results
%--------------------------------------------

%SWITCH CHANNEL TO LABEL?
grandavg.label     = cfg.channel;       % cell-array
grandavg.fsample   = varargin{1}.fsample;
grandavg.avg       = ResultGrandavg; 		% Nchan x Nsamples
grandavg.var       = ResultVar;		      % Nchan x Nsamples
grandavg.time      = ResultsTime;       % 1 x Nsamples

%KEEP INDIVIDUAL MEANS?
if strcmp(cfg.keepindividual, 'yes')
  grandavg.individual = avgmat;         % Nsubj x Nchan x Nsamples
  grandavg.dimord = 'subj_chan_time';
else
  grandavg.dimord = 'chan_time';
end

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
cfg.version.id = '$Id: timelockgrandaverage.m,v 1.21 2009/01/20 13:01:31 sashae Exp $';
% remember the configuration details of the input data
cfg.previous = [];
for i=1:length(varargin)
  try, cfg.previous{i} = varargin{i}.cfg; end
end
% remember the exact configuration details in the output
grandavg.cfg = cfg;

