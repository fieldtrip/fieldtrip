function [trl, event] = trialfun_general(cfg)

% TRIALFUN_GENERAL determines trials/segments in the data that are
% interesting for analysis, using the general event structure returned
% by read_event. This function is independent of the dataformat
%
% The trialdef structure can contain the following specifications
%   cfg.trialdef.eventtype  = 'string'
%   cfg.trialdef.eventvalue = number, string or list with numbers or strings
%   cfg.trialdef.prestim    = latency in seconds (optional)
%   cfg.trialdef.poststim   = latency in seconds (optional)
%
% If you want to read all data from a continous file in segments, you can specify
%    cfg.trialdef.triallength = duration in seconds (can be Inf)
%    cfg.trialdef.ntrials     = number of trials
%
% If you specify
%   cfg.trialdef.eventtype  = '?'
% a list with the events in your datafile will be displayed on screen.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

% Copyright (C) 2005-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% some events do not require the specification a type, pre or poststim period
% in that case it is more convenient not to have them, instead of making them empty
if isfield(cfg.trialdef, 'eventvalue')  && isempty(cfg.trialdef.eventvalue   ), cfg.trialdef = rmfield(cfg.trialdef, 'eventvalue' ); end
if isfield(cfg.trialdef, 'prestim')     && isempty(cfg.trialdef.prestim      ), cfg.trialdef = rmfield(cfg.trialdef, 'prestim'    ); end
if isfield(cfg.trialdef, 'poststim')    && isempty(cfg.trialdef.poststim     ), cfg.trialdef = rmfield(cfg.trialdef, 'poststim'   ); end
if isfield(cfg.trialdef, 'triallength') && isempty(cfg.trialdef.triallength  ), cfg.trialdef = rmfield(cfg.trialdef, 'triallength'); end
if isfield(cfg.trialdef, 'ntrials')     && isempty(cfg.trialdef.ntrials      ), cfg.trialdef = rmfield(cfg.trialdef, 'ntrials'    ); end

if isfield(cfg.trialdef, 'triallength')
  % reading all segments from a continuous fie is incompatible with any other option
  try, cfg.trialdef = rmfield(cfg.trialdef, 'eventvalue'); end
  try, cfg.trialdef = rmfield(cfg.trialdef, 'prestim'   ); end
  try, cfg.trialdef = rmfield(cfg.trialdef, 'poststim'  ); end
  if ~isfield(cfg.trialdef, 'ntrials')
    if isinf(cfg.trialdef.triallength)
      cfg.trialdef.ntrials = 1;
    else
      cfg.trialdef.ntrials = inf;
    end
  end
end

% default rejection parameter
if ~isfield(cfg, 'eventformat'),  cfg.eventformat  = []; end
if ~isfield(cfg, 'headerformat'), cfg.headerformat = []; end
if ~isfield(cfg, 'dataformat'),   cfg.dataformat   = []; end

% read the header, contains the sampling frequency 
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);

% read the events
if isfield(cfg, 'event')
  fprintf('using the events from the configuration structure\n');
  event = cfg.event;
else
  try
    fprintf('reading the events from ''%s''\n', cfg.headerfile);
    event = ft_read_event(cfg.headerfile, 'headerformat', cfg.headerformat, 'eventformat', cfg.eventformat, 'dataformat', cfg.dataformat);
  catch
    event = [];
  end
end

% for the following, the trials do not depend on the events in the data
if isfield(cfg.trialdef, 'triallength')
  if isinf(cfg.trialdef.triallength)
    % make one long trial with the complete continuous data in it
    trl = [1 hdr.nSamples*hdr.nTrials 0];
  elseif isinf(cfg.trialdef.ntrials)
    % cut the continous data into as many segments as possible
    nsamples = round(cfg.trialdef.triallength*hdr.Fs);
    trlbeg   = 1:nsamples:(hdr.nSamples*hdr.nTrials - nsamples + 1);
    trlend   = trlbeg + nsamples - 1;
    offset   = zeros(size(trlbeg));
    trl = [trlbeg(:) trlend(:) offset(:)];
  else
    % make the pre-specified number of trials
    nsamples = round(cfg.trialdef.triallength*hdr.Fs);
    trlbeg   = (0:(cfg.trialdef.ntrials-1))*nsamples + 1;
    trlend   = trlbeg + nsamples - 1;
    offset   = zeros(size(trlbeg));
    trl = [trlbeg(:) trlend(:) offset(:)];
  end
  return
end

sel = [];
trl = [];
val = [];
if strcmp(cfg.trialdef.eventtype, '?')
  % no trials should be added, show event information using subfunction and exit
  show_event(event);
  return
end

if strcmp(cfg.trialdef.eventtype, 'gui') || (isfield(cfg.trialdef, 'eventvalue') && length(cfg.trialdef.eventvalue)==1 && strcmp(cfg.trialdef.eventvalue, 'gui'))
  cfg.trialdef = select_event(event, cfg.trialdef);
  usegui = 1;
else
  usegui = 0;
end

if ~isfield(cfg.trialdef, 'eventvalue')
  cfg.trialdef.eventvalue = [];
elseif ischar(cfg.trialdef.eventvalue)
  % convert single string into cell-array, otherwise the intersection does not work as intended
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

% select all events of the specified type and with the specified value
if ~isempty(cfg.trialdef.eventtype)
  sel = ismember({event.type}, cfg.trialdef.eventtype);
else
  sel = true(size(event));
end

if ~isempty(cfg.trialdef.eventvalue)
  % this cannot be done robustly in a single line of code
  if ~iscell(cfg.trialdef.eventvalue)
    valchar    = ischar(cfg.trialdef.eventvalue); 
    valnumeric = isnumeric(cfg.trialdef.eventvalue); 
  else
    valchar    = ischar(cfg.trialdef.eventvalue{1});
    valnumeric = isnumeric(cfg.trialdef.eventvalue{1});
  end  
  for i=1:numel(event)
    if (ischar(event(i).value) && valchar) || (isnumeric(event(i).value) && valnumeric)
      sel(i) = sel(i) & ~isempty(intersect(event(i).value, cfg.trialdef.eventvalue));
    end
  end
end

% convert from boolean vector into a list of indices
sel = find(sel);  

if usegui
  % Checks whether offset and duration are defined for all the selected
  % events and/or prestim/poststim are defined in trialdef.
  if (any(cellfun('isempty', {event(sel).offset})) || ...
      any(cellfun('isempty', {event(sel).duration}))) && ...
      ~(isfield(cfg.trialdef, 'prestim') && isfield(cfg.trialdef, 'poststim'))
  
    % If at least some of offset/duration values and prestim/poststim
    % values are missing tries to ask the user for prestim/poststim
    answer = inputdlg({'Prestimulus latency (sec)','Poststimulus latency (sec)'}, 'Enter borders');
    if isempty(answer) || any(cellfun('isempty', answer))
      error('The information in the data and cfg is insufficient to define trials.');
    else
      cfg.trialdef.prestim=str2double(answer{1});
      cfg.trialdef.poststim=str2double(answer{2});
      if isnan(cfg.trialdef.prestim) || isnan(cfg.trialdef.poststim)
        error('Illegal input for trial borders');
      end
    end
  end % if specification is not complete
end % if usegui

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
  if ~isfield(cfg.trialdef, 'prestim')
    trloff = event(i).offset;
    trlbeg = event(i).sample;
  else
    % override the offset of the event
    trloff = round(-cfg.trialdef.prestim*hdr.Fs);
    % also shift the begin sample with the specified amount
    trlbeg = event(i).sample + trloff;
  end
  % determine the number of samples that has to be read (excluding the begin sample)
  if ~isfield(cfg.trialdef, 'poststim')
    trldur = max(event(i).duration - 1, 0);
  else
    % this will not work if prestim was not defined, the code will then crash
    trldur = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;
  end
  trlend = trlbeg + trldur;
  % add the beginsample, endsample and offset of this trial to the list
  % if all samples are in the dataset
  if trlbeg>0 && trlend<=hdr.nSamples.*hdr.nTrials,
    trl = [trl; [trlbeg trlend trloff]];
    if isnumeric(event(i).value), 
      val = [val; event(i).value];
    else
      val = [val; nan];
    end
  end
end

% append the vector with values
if ~isempty(val) && ~all(isnan(val))
  trl = [trl val];
end

if usegui
    % This complicated line just computes the trigger times in seconds and
    % converts them to a cell array of strings to use in the GUI
    eventstrings = cellfun(@num2str, mat2cell((trl(:, 1)- trl(:, 3))./hdr.Fs , ones(1, size(trl, 1))), 'UniformOutput', 0);

    % Let us start with handling at least the completely unsegmented case
    % semi-automatically. The more complicated cases are better left
    % to the user.
    if hdr.nTrials==1
        selected = find(trl(:,1)>0 & trl(:,2)<=hdr.nSamples);
    else
        selected = find(trl(:,1)>0);
    end

    indx = select_channel_list(eventstrings, selected , 'Select events');

    trl=trl(indx, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that shows event table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_event(event)
if isempty(event)
  fprintf('no events were found in the datafile\n');
  return
end
eventtype = unique({event.type});
Neventtype = length(eventtype);
if Neventtype==0
  fprintf('no events were found in the datafile\n');
else
  fprintf('the following events were found in the datafile\n');
  for i=1:Neventtype
    sel = find(strcmp(eventtype{i}, {event.type}));
    try
      eventvalue = unique({event(sel).value});            % cell-array with string value
      eventvalue = sprintf('''%s'' ', eventvalue{:});     % translate into a single string
    catch
      eventvalue = unique(cell2mat({event(sel).value}));  % array with numeric values or empty
      eventvalue = num2str(eventvalue);                   % translate into a single string
    end
    fprintf('event type: ''%s'' ', eventtype{i});
    fprintf('with event values: %s', eventvalue);
    fprintf('\n');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that allows the user to select an event using gui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trialdef=select_event(event, trialdef)
if isempty(event)
  fprintf('no events were found in the datafile\n');
  return
end
if strcmp(trialdef.eventtype, 'gui')
  eventtype = unique({event.type});
else
  eventtype ={trialdef.eventtype};
end
Neventtype = length(eventtype);
if Neventtype==0
  fprintf('no events were found in the datafile\n');
else
  % Two lists are built in parallel
  settings={}; % The list of actual values to be used later
  strsettings={}; % The list of strings to show in the GUI
  for i=1:Neventtype
    sel = find(strcmp(eventtype{i}, {event.type}));

    emptyval=find(cellfun('isempty', {event(sel).value}));

    if all(cellfun(@isnumeric, {event(sel).value}))
      [event(sel(emptyval)).value]=deal(Inf);
      eventvalue = unique([event(sel).value]);
    else
      if ~isempty(strmatch('Inf', {event(sel).value},'exact'))
        % It's a very unlikely scenario but ...
        warning('Event value''Inf'' cannot be handled by GUI selection. Mistakes are possible.')
      end
      [event(sel(emptyval)).value]=deal('Inf');
      eventvalue = unique({event(sel).value});
      if ~iscell(eventvalue)
        eventvalue={eventvalue};
      end
    end
    for j=1:length(eventvalue)
      if (ischar(eventvalue(j)) && ~strcmp(eventvalue(j), 'Inf')) || ...
          (isnumeric(eventvalue(j)) && eventvalue(j)~=Inf)
        settings=[settings; [eventtype(i), eventvalue(j)]];
      else
        settings=[settings; [eventtype(i), {[]}]];
      end

      if isa(eventvalue, 'numeric')
        strsettings=[strsettings; {['Type: ' eventtype{i} ' ; Value: ' num2str(eventvalue(j))]}];
      else
        strsettings=[strsettings; {['Type: ' eventtype{i} ' ; Value: ' eventvalue{j}]}];
      end
    end
  end
  if isempty(strsettings)
    fprintf('no events of the selected type were found in the datafile\n');
    return
  end

  [selection ok]= listdlg('ListString',strsettings, 'SelectionMode', 'single', 'Name', 'Select event', 'ListSize', [300 300]);

  if ok
    trialdef.eventtype=settings{selection,1};
    trialdef.eventvalue=settings{selection,2};
  end
end

