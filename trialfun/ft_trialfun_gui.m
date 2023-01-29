function [trl, event] = ft_trialfun_gui(cfg)

% FT_TRIALFUN_GUI reads events from the dataset, displays a graphical user interface
% dialog to select the event types and values of interest, and constructs a trial
% definition. This function should in general not be called directly, it will be
% called by FT_DEFINETRIAL.
%
% Use this function by calling
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset   = string with the filename
%   cfg.trialfun  = 'ft_trialfun_gui'
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL, FT_TRIALFUN_SHOW

% Copyright (C) 2005-2021, Robert Oostenveld, Vladimir Litvak
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

% set the defaults
cfg.trialdef.eventtype    = ft_getopt(cfg.trialdef, 'eventtype');
cfg.trialdef.eventvalue   = ft_getopt(cfg.trialdef, 'eventvalue');
cfg.trialdef.prestim      = ft_getopt(cfg.trialdef, 'prestim');
cfg.trialdef.poststim     = ft_getopt(cfg.trialdef, 'poststim');

% these options get passed to FT_READ_EVENT
cfg.trialdef.detectflank  = ft_getopt(cfg.trialdef, 'detectflank');
cfg.trialdef.trigshift    = ft_getopt(cfg.trialdef, 'trigshift');
cfg.trialdef.chanindx     = ft_getopt(cfg.trialdef, 'chanindx');
cfg.trialdef.threshold    = ft_getopt(cfg.trialdef, 'threshold');
cfg.trialdef.tolerance    = ft_getopt(cfg.trialdef, 'tolerance');

% specify the default file formats
cfg.eventformat   = ft_getopt(cfg, 'eventformat');
cfg.headerformat  = ft_getopt(cfg, 'headerformat');
cfg.dataformat    = ft_getopt(cfg, 'dataformat');

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
  event = ft_read_event(cfg.headerfile, 'headerformat', cfg.headerformat, 'eventformat', cfg.eventformat, 'dataformat', cfg.dataformat,  'detectflank', cfg.trialdef.detectflank, 'trigshift', cfg.trialdef.trigshift, 'chanindx', cfg.trialdef.chanindx, 'threshold', cfg.trialdef.threshold, 'tolerance', cfg.trialdef.tolerance);
end

% this function always returns an empty trial definition
trl = [];

if isempty(event)
  ft_info('no events were found\n');
  
else
  % make a pre-selection of event types
  if strcmp(cfg.trialdef.eventtype, 'gui') || isempty(cfg.trialdef.eventtype)
    selectedtypes = unique({event.type});
  elseif iscell(cfg.trialdef.eventtype)
    selectedtypes = cfg.trialdef.eventtype;
  elseif ischar(cfg.trialdef.eventtype)
    selectedtypes = {cfg.trialdef.eventtype};
  end
  
  % make a pre-selection of event values
  if strcmp(cfg.trialdef.eventvalue, 'gui') || isempty(cfg.trialdef.eventvalue)
    % simply concatenating all values fails if they are a mixture of numbers and strings
    selectedvalues = {event.value};
    sel = cellfun(@ischar, selectedvalues);
    uniquechar = unique(selectedvalues(sel));
    sel = cellfun(@isscalar, selectedvalues);
    uniquescalar = num2cell(unique(cell2mat(selectedvalues(sel))));
    % determine the combination of both unique sets
    selectedvalues = [uniquescalar uniquechar];
  elseif iscell(cfg.trialdef.eventvalue)
    selectedvalues = cfg.trialdef.eventvalue;
  elseif isnumeric(cfg.trialdef.eventvalue)
    selectedvalues = num2cell(cfg.trialdef.eventvalue);
  end
  
  % Two lists are built, the list of real values to be used later and the list of strings to show in the GUI
  listreal = {};
  listshow = {};
  for i=1:length(event)
    if ~ismatch(event(i).type, selectedtypes)
      continue
    elseif ~isempty(event(i).value) && ~ismatch(event(i).value, selectedvalues)
      continue
    end
    
    listreal{end+1,1} = event(i).type;
    listreal{end,  2} = event(i).value;
    
    if isempty(event(i).value)
      listshow{end+1} = ['Type: ' event(i).type ' ; Value: []'];
    elseif isnumeric(event(i).value)
      listshow{end+1} = ['Type: ' event(i).type ' ; Value: ' num2str(event(i).value)];
    else
      listshow{end+1} = ['Type: ' event(i).type ' ; Value: ' event(i).value];
    end
  end % for each event
  
  % only keep the unique combinations
  [listshow, index] = unique(listshow);
  listreal = listreal(index,:);
  
  % show the selection dialog with the event types and values
  [selection, ok] = listdlg('ListString', listshow, 'SelectionMode', 'multiple', 'Name', 'Select event', 'ListSize', [300 300]);
  
  if ok
    cfg.trialdef.eventtype  = listreal(selection,1);
    cfg.trialdef.eventvalue = listreal(selection,2);
    % empty values are not interesting
    sel = cellfun(@isempty, cfg.trialdef.eventvalue);
    cfg.trialdef.eventvalue(sel) = [];
  else
    trl = [];
    return
  end
  
  % make the initial selection
  sel = true(size(event));
  
  % select all events of the specified type
  if ~isempty(cfg.trialdef.eventtype)
    for i=1:numel(event)
      sel(i) = sel(i) && ismatch(event(i).type, cfg.trialdef.eventtype);
    end
  end
  
  % select all events with the specified value
  if ~isempty(cfg.trialdef.eventvalue)
    for i=1:numel(event)
      sel(i) = sel(i) && ismatch(event(i).value, cfg.trialdef.eventvalue);
    end
  end
  
  % convert from boolean vector into a list of indices, this must be a row vector
  sel = find(sel(:)');
  
  % Check whether offset/duration and/or prestim/poststim are specified
  if (any(cellfun('isempty', {event(sel).offset})) || ...
      any(cellfun('isempty', {event(sel).duration}))) && ...
      (isempty(cfg.trialdef.prestim) || isempty(cfg.trialdef.poststim))
    
    % If at least some of offset/duration values and prestim/poststim
    % values are missing we will ask the user for prestim/poststim
    answer = inputdlg({'Prestimulus latency (sec)','Poststimulus latency (sec)'}, 'Enter borders');
    if isempty(answer) || any(cellfun('isempty', answer))
      ft_error('The information in the data and cfg is insufficient to define trials.');
    else
      cfg.trialdef.prestim  = str2double(answer{1});
      cfg.trialdef.poststim = str2double(answer{2});
      if isnan(cfg.trialdef.prestim) || isnan(cfg.trialdef.poststim)
        ft_error('Illegal input for trial borders');
      end
    end
  end % if the specification is not complete
  
  for i=sel
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
    elseif ischar(event(i).value) && numel(event(i).value)>1 && (event(i).value(1)=='S'|| event(i).value(1)=='R')
      % on brainvision these are called 'S  1' for stimuli or 'R  1' for responses
      trlval = str2double(event(i).value(2:end));
    else
      trlval = nan;
    end
    
    % add the trial only if all samples are in the dataset
    if trlbeg>0 && trlend<=hdr.nSamples*hdr.nTrials
      thistrl = [trlbeg trlend trloff trlval];
      trl = cat(1, trl, thistrl);
    end
  end
  
  if ~isempty(trl) && all(isnan(trl(:,4)))
    % the values are not informative, remove them
    trl = trl(:,1:3);
  end
  
  % This complicated line just computes the trial onset times in seconds
  % and converts them to a cell-array of strings to use in the GUI
  triggertime = cellfun(@num2str, mat2cell((trl(:, 1)-trl(:, 3))./hdr.Fs, ones(1, size(trl, 1))), 'UniformOutput', 0);
  
  % Let us start with handling at least the completely unsegmented case
  % semi-automatically. The more complicated cases are better left to the user.
  if hdr.nTrials==1
    selected = find(trl(:,1)>0 & trl(:,2)<=hdr.nSamples);
  else
    selected = find(trl(:,1)>0);
  end
  
  indx = select_channel_list(triggertime, selected, 'Select event latency');
  
  trl = trl(indx, :);
  
end % if isempty(event)
