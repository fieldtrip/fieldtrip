function [data] = redefinetrial(cfg, data)

% REDEFINETRIAL allows you to adjust the time axis of your data, i.e. to
% change from stimulus-locked to response-locked. Furthermore, it allows
% you to select a time window of interest, or to resegment your long trials
% into shorter fragments.
%
% Use as
%   data = redefinetrial(cfg, data)
% where the input data should correspond to the output of PREPROCESSING and
% the configuration should be specified as explained below. Note that some
% options are mutually exclusive, and require two calls to this function to
% avoid confucion about the order in which they are applied.
%
% For selecting a subset of trials you can specify
%   cfg.trials    = 'all' or a selection given as a 1xN vector (default = 'all')
%
% For selecting trials with a minimum length you can specify
%   cfg.minlength = length in seconds, can be 'maxperlen' (default = [])
%
% For realiging the time axes of all trials to a new reference time
% point (i.e. change the definition for t=0) you can use the following
% configuration option
%   cfg.offset    = single number or Nx1 vector, expressed in samples relative to current t=0
%
% For selecting a specific subsection of (i.e. cut out a time window
% of interest) you can select a time window in seconds that is common
% in all trials
%   cfg.toilim    = [tmin tmax] to specify a latency window in seconds
%
% Alternatively you can specify the begin and end sample in each trial
%   cfg.begsample = single number or Nx1 vector, expressed in samples relative to the start of the input trial
%   cfg.endsample = single number or Nx1 vector, expressed in samples relative to the start of the input trial
%
% Alternatively you can specify a new trial definition, expressed in
% samples relative to the original recording
%   cfg.trl       = Nx3 matrix with the trial definition, see DEFINETRIAL
%
% See also DEFINETRIAL, RECODEEVENT, PREPROCESSING

% Copyright (C) 2006-2008, Robert Oostenveld
%
% $Log: redefinetrial.m,v $
% Revision 1.25  2009/08/13 20:47:02  jansch
% changed typo in key-value pair 'hdr',hdr. this should read 'header',hdr and
% speeds up the function considerably
%
% Revision 1.24  2009/01/14 15:11:22  roboos
% corrected documentation for cfg.begsample/endsample
% fixed bug in output data.cfg.trl (one sample too long at the end) for input cfg.toilim and input cfg.begsample/endsample
% cleaned up the code for cfg.trl (use default variable names)
%
% Revision 1.23  2008/11/10 10:48:55  jansch
% removed typo
%
% Revision 1.22  2008/11/10 10:47:16  jansch
% added 1 to recomputed time-axes if cfg.trl. before, the time-axis was one
% sample too short
%
% Revision 1.21  2008/11/10 10:01:58  jansch
% ensured correct trial-definition in output data.cfg.trl if cfg.trl
%
% Revision 1.20  2008/11/10 09:41:21  jansch
% added sampling-frequency to output if cfg.trl is specified
%
% Revision 1.19  2008/10/07 16:08:11  estmee
% Changed the way fetch_data is called.
%
% Revision 1.18  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.17  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.16  2008/07/30 07:43:40  roboos
% determine max sample number correctly, also if trials are ordered differently, i.e. trl(end,:) does not have to be the last one in the recording
%
% Revision 1.15  2008/06/25 06:37:57  roboos
% change in whitespace
%
% Revision 1.14  2008/05/06 14:03:10  sashae
% change in trial selection, cfg.trials can be a logical
%
% Revision 1.13  2008/04/21 14:27:12  jansch
% changed trlnew into trl to be able to output the new trl-matrix correctly
% and elegantly into the output-structure
%
% Revision 1.12  2008/04/15 16:11:54  estmee
% updated documentation
%
% Revision 1.11  2008/04/15 16:07:24  estmee
% inserted code that creates new trials based on cfg.trl
%
% Revision 1.10  2008/04/14 14:59:31  roboos
% added test for correct cfg.trl
%
% Revision 1.9  2008/04/08 20:21:32  roboos
% updated documentation, made some room in the code to plug in the new functionality from Esther
%
% Revision 1.8  2007/12/18 16:21:02  sashae
% added option for trial selection, updated documentation and fixed some old code,
% now uses findcfg for finding the trl
%
% Revision 1.7  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.6  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.5  2007/02/07 12:09:56  roboos
% fixed some small bugs, thanks to Nienke
%
% Revision 1.4  2007/01/04 12:18:05  roboos
% changed a print statement
%
% Revision 1.3  2007/01/02 09:48:06  roboos
% fixed some small typos, prevent the selection of data that is completely outside teh toilim, implemented support for cfg.minlength
%
% Revision 1.2  2006/10/24 06:48:45  roboos
% small changes and bugfixes to make it work, updated documentation
%
% Revision 1.1  2006/10/19 16:06:51  roboos
% new implementation
%

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'offset'),    cfg.offset = [];    end
if ~isfield(cfg, 'toilim'),    cfg.toilim = [];    end
if ~isfield(cfg, 'begsample'), cfg.begsample = []; end
if ~isfield(cfg, 'endsample'), cfg.endsample = []; end
if ~isfield(cfg, 'minlength'), cfg.minlength = []; end
if ~isfield(cfg, 'trials'),    cfg.trials = 'all'; end
if ~isfield(cfg, 'feedback'),  cfg.feedback = 'yes'; end
if ~isfield(cfg, 'trl'),       cfg.trl =  [];      end

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', cfg.feedback);
fb   = strcmp(cfg.feedback, 'yes');

% trl is not specified in the function call, but the data is given ->
% try to locate the trial definition (trl) in the nested configuration
if isfield(data,'cfg')
  trl = findcfg(data.cfg, 'trl');
  if length(data.trial)~=size(trl,1) || length(data.time)~=size(trl,1)
    error('the trial definition is inconsistent with the data');
  end
else
  trl = [];
end
trlold = trl;

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
  if fb, fprintf('selecting %d trials\n', length(cfg.trials)); end
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
  trl         = trl(cfg.trials,:);
  if length(cfg.offset)>1 && length(cfg.offset)~=length(cfg.trials)
    cfg.offset=cfg.offset(cfg.trials);
  end
  if length(cfg.begsample)>1 && length(cfg.begsample)~=length(cfg.trials)
    cfg.begsample=cfg.begsample(cfg.trials);
  end
  if length(cfg.endsample)>1 && length(cfg.endsample)~=length(cfg.trials)
    cfg.endsample=cfg.endsample(cfg.trials);
  end
end
Ntrial = size(trl,1);

% check the input arguments, only one method for processing is allowed
numoptions = ~isempty(cfg.toilim) + ~isempty(cfg.offset) + (~isempty(cfg.begsample) || ~isempty(cfg.endsample)) + ~isempty(cfg.trl);
if numoptions>1
  error('you should specify only one of the options for redefining the data segments');
end
if numoptions==0 && isempty(cfg.minlength) && strcmp(cfg.trials, 'all')
  error('you should specify at least one configuration option');
end

% start processing
if ~isempty(cfg.toilim)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % select a latency window from each trial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  begsample = zeros(Ntrial,1);
  endsample = zeros(Ntrial,1);
  skiptrial = zeros(Ntrial,1);
  for i=1:Ntrial
    if cfg.toilim(1)>data.time{i}(end) || cfg.toilim(2)<data.time{i}(1)
      begsample(i) = nan;
      endsample(i) = nan;
      skiptrial(i) = 1;
    else
      begsample(i) = nearest(data.time{i}, cfg.toilim(1));
      endsample(i) = nearest(data.time{i}, cfg.toilim(2));
      data.trial{i} = data.trial{i}(:, begsample(i):endsample(i));
      data.time{i}  = data.time{i} (   begsample(i):endsample(i));
    end
  end
  % also correct the trial definition
  if ~isempty(trl)
    trl(:,1) = trl(:,1) + begsample - 1;
    trl(:,2) = trl(:,1) + endsample - begsample;
    trl(:,3) = trl(:,3) + begsample - 1;
  end

  % remove trials that are completely empty
  trl = trl(~skiptrial,:);
  data.time  = data.time(~skiptrial);
  data.trial = data.trial(~skiptrial);
  if fb, fprintf('removing %d trials in which no data was selected\n', sum(skiptrial)); end

elseif ~isempty(cfg.offset)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % shift the time axis from each trial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  offset = cfg.offset(:);
  if length(cfg.offset)==1
    offset = repmat(offset, Ntrial, 1);
  end
  for i=1:Ntrial
    data.time{i} = data.time{i} + offset(i)/data.fsample;
  end

  % also correct the trial definition
  if ~isempty(trl)
    trl(:,3) = trl(:,3) + offset;
  end

elseif ~isempty(cfg.begsample) || ~isempty(cfg.endsample)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % select a latency window from each trial based on begin and/or end sample
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  begsample = cfg.begsample(:);
  endsample = cfg.endsample(:);
  if length(begsample)==1
    begsample = repmat(begsample, Ntrial, 1);
  end
  if length(endsample)==1
    endsample = repmat(endsample, Ntrial, 1);
  end
  for i=1:Ntrial
    data.trial{i} = data.trial{i}(:, begsample(i):endsample(i));
    data.time{i}  = data.time{i} (   begsample(i):endsample(i));
  end

  % also correct the trial definition
  if ~isempty(trl)
    trl(:,1) = trl(:,1) + begsample - 1;
    trl(:,2) = trl(:,1) + endsample - begsample;
    trl(:,3) = trl(:,3) + begsample - 1;
  end

elseif ~isempty(cfg.trl)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % select new trials from the existing data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dataold = data;   % make a copy of the old data
  clear data        % this line is very important, we want to completely reconstruct the data from the old data!

  % make header
  hdr = fetch_header(dataold);

  % make new data structure
  trl = cfg.trl;
  for iTrl=1:length(trl(:,1))
    begsample = trl(iTrl,1);
    endsample = trl(iTrl,2);
    offset    = trl(iTrl,3);
    trllength        = endsample - begsample + 1;
    data.trial{iTrl} = fetch_data(dataold, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', 1:hdr.nChans, 'docheck', 0);
    data.time{iTrl}  = offset2time(offset, dataold.fsample, trllength);
  end
  data.hdr       = hdr;
  data.label     = dataold.label;
  data.fsample   = dataold.fsample;
  if isfield(dataold, 'grad')
    data.grad      = dataold.grad;
  end
  if isfield(dataold, 'elec')
    data.elec      = dataold.elec;
  end
end % processing the realignment or data selection

if ~isempty(cfg.minlength)
  Ntrial    = length(data.trial);
  trllength = zeros(Ntrial, 1);
  % determine the length of each trial
  for i=1:Ntrial
    trllength(i) = data.time{i}(end) - data.time{i}(1);
  end
  if ischar(cfg.minlength) && strcmp(cfg.minlength, 'maxperlen')
    minlength = max(trllength);
  else
    minlength = cfg.minlength;
  end
  % remove trials that are too short
  skiptrial = (trllength<minlength);
  trl = trl(~skiptrial,:);
  data.time  = data.time(~skiptrial);
  data.trial = data.trial(~skiptrial);
  if fb, fprintf('removing %d trials that are too short\n', sum(skiptrial)); end
end

% remember the previous and the up-to-date trial definitions in the configuration
cfg.trl    = trl;
cfg.trlold = trlold;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: redefinetrial.m,v 1.25 2009/08/13 20:47:02 jansch Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg;    end
try, cfg.previous = dataold.cfg; end % in case of ~isempty(cfg.trl)
% remember the exact configuration details in the output
data.cfg = cfg;

