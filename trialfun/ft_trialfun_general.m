function [trl, event] = ft_trialfun_general(cfg)

% FT_TRIALFUN_GENERAL reads events from the dataset using FT_READ_EVENT and
% constructs a trial definition. This function should in general not be called
% directly, it will be called by FT_DEFINETRIAL.
%
% Use this function by calling
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset   = string with the filename
%   cfg.trialdef  = structure with the details of trial definition, see below
%   cfg.trialfun  = 'ft_trialfun_general'
%
% The cfg.trialdef structure can contain the following specifications
%   cfg.trialdef.eventtype  = string, or cell-array with strings
%   cfg.trialdef.eventvalue = number, string, or list with numbers or strings
%   cfg.trialdef.prestim    = number, latency in seconds (optional)
%   cfg.trialdef.poststim   = number, latency in seconds (optional)
%
% You can specify these options that are passed to FT_READ_EVENT for trigger detection
%   cfg.trialdef.detectflank  = string, can be 'up', 'updiff', 'down', 'downdiff', 'both', 'any', 'biton', 'bitoff'
%   cfg.trialdef.trigshift    = integer, number of samples to shift from flank to detect trigger value
%   cfg.trialdef.chanindx     = list with channel numbers for the trigger detection, specify -1 in case you don't want to detect triggers
%   cfg.trialdef.threshold    = threshold for analog trigger channels
%   cfg.trialdef.tolerance    = tolerance in samples when merging analogue trigger channels, only for Neuromag
%
% If you want to read all data from a continuous file in segments, you can specify
%    cfg.trialdef.length      = duration of the segments in seconds (can be Inf)
%    cfg.trialdef.ntrials     = number of trials (optional, can be 1)
%    cfg.trialdef.overlap     = single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GUI, FT_TRIALFUN_SHOW

% Copyright (C) 2005-2021, Robert Oostenveld
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

% most defaults are in trialdef
cfg.trialdef = ft_getopt(cfg, 'trialdef', struct());

% check if the input cfg is valid for this function
cfg.trialdef = ft_checkconfig(cfg.trialdef, 'renamed', {'triallength', 'length'});
cfg.trialdef = ft_checkconfig(cfg.trialdef, 'renamedval', {'ntrials', inf, []});

% set the defaults
cfg.trialdef.eventtype    = ft_getopt(cfg.trialdef, 'eventtype');
cfg.trialdef.eventvalue   = ft_getopt(cfg.trialdef, 'eventvalue');
cfg.trialdef.prestim      = ft_getopt(cfg.trialdef, 'prestim');
cfg.trialdef.poststim     = ft_getopt(cfg.trialdef, 'poststim');

% these options are similar to those in FT_REDEFINETRIALS
cfg.trialdef.length       = ft_getopt(cfg.trialdef, 'length');
cfg.trialdef.overlap      = ft_getopt(cfg.trialdef, 'overlap', 0); % between 0 and 1
cfg.trialdef.ntrials      = ft_getopt(cfg.trialdef, 'ntrials');

% these options get passed to FT_READ_EVENT
cfg.trialdef.detectflank  = ft_getopt(cfg.trialdef, 'detectflank');
cfg.trialdef.trigshift    = ft_getopt(cfg.trialdef, 'trigshift');
cfg.trialdef.chanindx     = ft_getopt(cfg.trialdef, 'chanindx');
cfg.trialdef.threshold    = ft_getopt(cfg.trialdef, 'threshold');
cfg.trialdef.tolerance    = ft_getopt(cfg.trialdef, 'tolerance');
cfg.trialdef.combinebinary = ft_getopt(cfg.trialdef, 'combinebinary');

% specify the default file formats
cfg.eventformat   = ft_getopt(cfg, 'eventformat');
cfg.headerformat  = ft_getopt(cfg, 'headerformat');
cfg.dataformat    = ft_getopt(cfg, 'dataformat');
cfg.representation = ft_getopt(cfg, 'representation');

% get the header, this is among others for the sampling frequency
if isfield(cfg, 'hdr')
  ft_info('using the header from the configuration structure\n');
  hdr = cfg.hdr;
else
  % read the header, contains the sampling frequency
  ft_info('reading the header from ''%s''\n', cfg.headerfile);
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
end

% get the events
if isfield(cfg, 'event')
  ft_info('using the events from the configuration structure\n');
  event = cfg.event;
else
  ft_info('reading the events from ''%s''\n', cfg.headerfile);
  event = ft_read_event(cfg.headerfile, 'headerformat', cfg.headerformat, 'eventformat', cfg.eventformat, 'dataformat', cfg.dataformat,  'detectflank', cfg.trialdef.detectflank, 'trigshift', cfg.trialdef.trigshift, 'chanindx', cfg.trialdef.chanindx, 'threshold', cfg.trialdef.threshold, 'tolerance', cfg.trialdef.tolerance, 'combinebinary', cfg.trialdef.combinebinary);
end

if ~isempty(cfg.trialdef.length) && ~isinf(cfg.trialdef.length)
  % make as many trials as possible with the specified length and offset
  begsample   = 1;
  endsample   = round(hdr.nSamples*hdr.nTrials);
  offset      = 0;
  nsmp        = round(cfg.trialdef.length*hdr.Fs);
  nshift      = round((1-cfg.trialdef.overlap)*nsmp);
  alltrl      = (begsample:nshift:(endsample+1-nsmp))';
  alltrl(:,2) = alltrl(:,1) + nsmp - 1;
  alltrl(:,3) = alltrl(:,1) + offset - alltrl(1,1);
  % trim to the requested number of trials
  if ~isempty(cfg.trialdef.ntrials)
    trl = alltrl(1:cfg.trialdef.ntrials,:);
  else
    trl = alltrl;
  end
  return
  
elseif isscalar(cfg.trialdef.ntrials) || isequal(cfg.trialdef.length, inf)
  % construct a single trial
  if isscalar(cfg.trialdef.ntrials) && cfg.trialdef.ntrials~=1
    ft_error('this is only supported for a single trial');
  end
  begsample = 1;
  endsample = round(hdr.nSamples*hdr.nTrials);
  offset    = 0;
  trl       = [begsample endsample offset];
  return
  
else
  % select events on basis of event types and values
  sel = true(1, length(event)); % this should be a row vector
  
  % select all events of the specified type
  if isfield(cfg.trialdef, 'eventtype') && ~isempty(cfg.trialdef.eventtype)
    for i=1:numel(event)
      sel(i) = sel(i) && ismatch(event(i).type, cfg.trialdef.eventtype);
    end
  elseif isempty(cfg.trialdef.eventtype)
    % search for trial events
    for i=1:numel(event)
      sel(i) = sel(i) && ismatch(event(i).type, 'trial');
    end
  end
  
  % select all events with the specified value
  if isfield(cfg.trialdef, 'eventvalue') && ~isempty(cfg.trialdef.eventvalue)
    for i=1:numel(event)
      sel(i) = sel(i) && ismatch(event(i).value, cfg.trialdef.eventvalue);
    end
  end
  
  % convert from boolean vector into a list of indices
  sel = find(sel);
  
  % start with an empty list
  trl = [];
  
  for i=sel
    % catch empty fields in the event table and interpret them meaningfully
    if isempty(event(i).offset)
      % time axis has no offset relative to the event
      event(i).offset = 0;
    end
    if isempty(event(i).duration)
      % the event does not specify a duration
      event(i).duration = 0;
    end
    
    % determine where the trial starts with respect to the event
    if isempty(cfg.trialdef.prestim)
      trloff = event(i).offset;
      trlbeg = event(i).sample;
    else
      % override the offset of the event
      trloff = round(-cfg.trialdef.prestim*hdr.Fs);
      % also shift the begin sample with the specified amount
      trlbeg = event(i).sample + trloff;
    end
    % determine the number of samples that has to be read (excluding the begin sample)
    if isempty(cfg.trialdef.poststim)
      trldur = max(event(i).duration - 1, 0);
    else
      % this will not work if prestim was not defined, the code will then crash
      trldur = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;
    end
    trlend = trlbeg + trldur;
    
    if isnumeric(event(i).value) && ~isempty(event(i).value)
      trlval = event(i).value;
    elseif ischar(event(i).value) && ~isempty(regexp(event(i).value, '^[SR]+[\s]*+[0-9]{1,3}$'))
      % This looks like Brainvision event markers. For backward compatibility, convert
      % the strings into the numerals following the 'S' or 'R', unless the user has specified
      % the cfg.representation to be a table
      if ~isequal(cfg.representation, 'table')
        ft_warning('Brainvision markers are converted to numeric representation, if you want tabular output please specify cfg.representation=''table''');
        trlval = str2double(event(i).value(2:end));
      else
        trlval = event(i).value;
      end
    elseif ischar(event(i).value) && ~isequal(cfg.representation, 'numeric')
      trlval = event(i).value;
    else
      % the following depends on cfg.representation
      if isequal(cfg.representation, 'numeric') || isempty(cfg.representation)
        trlval = nan;
      else
        trlval = event(i).value;
      end
    end
    
    % add the trial only if all samples are in the dataset
    if trlbeg>0 && trlend<=hdr.nSamples*hdr.nTrials
      if isnumeric(trlval)
        % create a numeric array
        thistrl = [trlbeg trlend trloff trlval];
      else
        thistrl = cell2table({trlbeg trlend trloff trlval});
      end
      trl = cat(1, trl, thistrl);
    end
  end
  
  if ~isempty(trl) && ~istable(trl) && all(isnan(trl(:,4)))
    % the values are not informative, remove them
    trl = trl(:,1:3);
  elseif ~isempty(trl) && istable(trl)
    % add names to the columns of the table
    trl.Properties.VariableNames = {'begsample', 'endsample', 'offset', 'eventvalue'};
  end
  
  
end
