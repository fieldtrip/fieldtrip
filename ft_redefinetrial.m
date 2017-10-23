function [data] = ft_redefinetrial(cfg, data)

% FT_REDEFINETRIAL allows you to adjust the time axis of your data, i.e. to
% change from stimulus-locked to response-locked. Furthermore, it allows
% you to select a time window of interest, or to resegment your long trials
% into shorter fragments.
%
% Use as
%   data = ft_redefinetrial(cfg, data)
% where the input data should correspond to the output of FT_PREPROCESSING and
% the configuration should be specified as explained below. Note that some
% options are mutually exclusive, and require two calls to this function to
% avoid confusion about the order in which they are applied.
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
%   cfg.trl       = Nx3 matrix with the trial definition, see FT_DEFINETRIAL
%
% Alternatively you can specify the data to be cut into (non-)overlapping
% segments, starting from the beginning of each trial. This may lead to loss
% of data at the end of the trials
%   cfg.length    = single number (in unit of time, typically seconds) of the required snippets
%   cfg.overlap   = single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_DEFINETRIAL, FT_RECODEEVENT, FT_PREPROCESSING

% Copyright (C) 2006-2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% ft_checkdata is done further down

% set the defaults
cfg.offset       = ft_getopt(cfg, 'offset',    []);
cfg.toilim       = ft_getopt(cfg, 'toilim',    []);
cfg.begsample    = ft_getopt(cfg, 'begsample', []);
cfg.endsample    = ft_getopt(cfg, 'endsample', []);
cfg.minlength    = ft_getopt(cfg, 'minlength', []);
cfg.trials       = ft_getopt(cfg, 'trials',    'all', 1);
cfg.feedback     = ft_getopt(cfg, 'feedback',  'yes');
cfg.trl          = ft_getopt(cfg, 'trl',       []);
cfg.length       = ft_getopt(cfg, 'length',    []);
cfg.overlap      = ft_getopt(cfg, 'overlap',   0);

% store original datatype
dtype = ft_datatype(data);

% deal with the special case of timelock rpt_chan_time with 1 trial
oneRptTimelock = (strcmp(dtype, 'timelock') &&...
  strcmp(data.dimord, 'rpt_chan_time') &&...
  size(data.trial, 1) == 1);


% check if the input data is valid for this function, this will convert it to raw if needed
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', cfg.feedback);
fb   = istrue(cfg.feedback);

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  if fb, fprintf('selecting %d trials\n', length(cfg.trials)); end

  % select trials of interest
  tmpcfg = keepfields(cfg, {'trials', 'showcallinfo'});
  data   = ft_selectdata(tmpcfg, data);
  % restore the provenance information
  [cfg, data] = rollback_provenance(cfg, data);

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
Ntrial = numel(data.trial);

% check the input arguments, only one method for processing is allowed
numoptions = ~isempty(cfg.toilim) + ~isempty(cfg.offset) + (~isempty(cfg.begsample) || ~isempty(cfg.endsample)) + ~isempty(cfg.trl) + ~isempty(cfg.length);
if numoptions>1
  ft_error('you should specify only one of the options for redefining the data segments');
end
if numoptions==0 && isempty(cfg.minlength) && strcmp(cfg.trials, 'all')
  ft_error('you should specify at least one configuration option');
end

% start processing
if ~isempty(cfg.toilim)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % select a latency window from each trial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  begsample = zeros(Ntrial,1);
  endsample = zeros(Ntrial,1);
  skiptrial = false(Ntrial,1);
  for i=1:Ntrial
    if cfg.toilim(1)>data.time{i}(end) || cfg.toilim(2)<data.time{i}(1)
      begsample(i) = nan;
      endsample(i) = nan;
      skiptrial(i) = true;
    else
      begsample(i) = nearest(data.time{i}, cfg.toilim(1));
      endsample(i) = nearest(data.time{i}, cfg.toilim(2));
      data.trial{i} = data.trial{i}(:, begsample(i):endsample(i));
      data.time{i}  = data.time{i} (   begsample(i):endsample(i));
    end
  end

  % also correct the sample information
  if isfield(data, 'sampleinfo')
    data.sampleinfo(:, 1) = data.sampleinfo(:, 1) + begsample - 1;
    data.sampleinfo(:, 2) = data.sampleinfo(:, 1) + endsample - begsample;
  end

  data.time     = data.time(~skiptrial);
  data.trial    = data.trial(~skiptrial);
  if isfield(data, 'sampleinfo'),  data.sampleinfo  = data.sampleinfo(~skiptrial, :); end
  if isfield(data, 'trialinfo'),   data.trialinfo   = data.trialinfo(~skiptrial, :);  end
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

  % also correct the sampleinfo
  if isfield(data, 'sampleinfo')
    data.sampleinfo(:, 1) = data.sampleinfo(:, 1) + begsample - 1;
    data.sampleinfo(:, 2) = data.sampleinfo(:, 1) + endsample - begsample;
  end

elseif ~isempty(cfg.trl)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % select new trials from the existing data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % ensure that sampleinfo is present, otherwise ft_fetch_data will crash
  data = ft_checkdata(data, 'hassampleinfo', 'yes');

  dataold = data;   % make a copy of the old data
  clear data        % this line is very important, we want to completely reconstruct the data from the old data!

  % make header
  hdr = ft_fetch_header(dataold);

  trl = cfg.trl;

  % start with a completely new data structure
  data          = [];
  data.hdr      = hdr;
  data.trial    = cell(1,size(trl,1));
  data.time     = cell(1,size(trl,1));
  data          = copyfields(dataold, data, {'fsample' 'label' 'topo' 'topolabel' 'unmixing' 'mixing' 'grad' 'elec' 'opto'}); % account for all potential fields to be copied over

  if isfield(dataold,'trialinfo')
    ft_warning('Original data has trialinfo, using user specified trialinfo instead');
  end
  
  for iTrl=1:length(trl(:,1))
    begsample = trl(iTrl,1);
    endsample = trl(iTrl,2);
    offset    = trl(iTrl,3);
    trllength = endsample - begsample + 1;

    data.trial{iTrl} = ft_fetch_data(dataold, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', 1:hdr.nChans, 'skipcheckdata', 1);
    data.time{iTrl}  = offset2time(offset, dataold.fsample, trllength);

    % ensure correct handling of trialinfo.
    % original trial
    iTrlorig  =  find(begsample <= dataold.sampleinfo(:,2) & endsample >= dataold.sampleinfo(:,1)); % Determines which old trials are present in new trials
    
    if size(cfg.trl,2)>3 %In case user specified a trialinfo
      data.trialinfo(iTrl,:) = cfg.trl(iTrl,4:end);
    elseif isfield(dataold,'trialinfo') % If old data has trialinfo
      if (numel(iTrlorig) == 1 ...  % only 1 old trial to copy trialinfo from, or
          || size(unique(dataold.trialinfo(iTrlorig,:),'rows'),1)) ... % all old trialinfo rows are identical
          && ~any(diff(dataold.sampleinfo(:,1))<=0) % and the trials are consecutive segments
        data.trialinfo(iTrl,:) = dataold.trialinfo(iTrlorig(1),:);
      else
        ft_error('Old trialinfo cannot be combined into new trialinfo, please specify trialinfo in cfg.trl(:,4)');
      end
    end
  end %for iTrl

  % adjust the sampleinfo in the output
  if isfield(dataold, 'sampleinfo')
    % adjust the trial definition
    data.sampleinfo  = trl(:, 1:2);
  end

elseif ~isempty(cfg.length)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % cut the existing trials into segments of the specified length
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  data = ft_checkdata(data, 'hassampleinfo', 'yes');

  % create dummy trl-matrix and recursively call ft_redefinetrial
  nsmp    = round(cfg.length*data.fsample);
  nshift  = round((1-cfg.overlap)*nsmp);

  newtrl = zeros(0,3);
  for k = 1:numel(data.trial)
    offset = time2offset(data.time{k}, data.fsample);
    tmp1   = [data.sampleinfo(k,:) offset];
    tmp2   = (tmp1(1):nshift:(tmp1(2)+1-nsmp))';
    if ~isempty(tmp2)
      tmp2(:,2) = tmp2 + nsmp - 1;
      tmp2(:,3) = tmp2(:,1) + offset - tmp2(1,1);
      newtrl = [newtrl; tmp2];
    end
  end

  tmpcfg = keepfields(cfg, {'showcallinfo'});
  tmpcfg.trl = newtrl;
  data   = ft_redefinetrial(tmpcfg, data);
  % restore the provenance information
  [cfg, data] = rollback_provenance(cfg, data);

end % processing the realignment or data selection

if ~isempty(cfg.minlength)
  Ntrial    = length(data.trial);
  trllength = zeros(Ntrial, 1);
  % determine the length of each trial
  for i=1:Ntrial
    trllength(i) = size(data.trial{i},2) * 1/data.fsample; % this the the DURATION of the selected samples
  end
  if ischar(cfg.minlength) && strcmp(cfg.minlength, 'maxperlen')
    minlength = max(trllength);
  else
    minlength = cfg.minlength;
  end
  % remove trials that are too short
  skiptrial = (trllength<minlength);
  %if ~isempty(trl), trl = trl(~skiptrial,:); end
  data.time  = data.time(~skiptrial);
  data.trial = data.trial(~skiptrial);
  if isfield(data, 'sampleinfo'), data.sampleinfo  = data.sampleinfo(~skiptrial, :); end
  if isfield(data, 'trialinfo'),  data.trialinfo   =  data.trialinfo(~skiptrial, :); end
  if fb, fprintf('removing %d trials that are too short\n', sum(skiptrial));         end
end

% convert back to input type if necessary
switch dtype
  case 'timelock'
    data = ft_checkdata(data, 'datatype', 'timelock');
    if oneRptTimelock
      % deal with the special case of rpt_chan_time timelock data with one
      % repetition
      data.trial = reshape(data.avg, [1 size(data.avg)]);
      data.dimord = 'rpt_chan_time';
    end
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
if ~isempty(cfg.trl)
  % the input data has been renamed to dataold
  ft_postamble previous dataold
else
  ft_postamble previous data
end
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
