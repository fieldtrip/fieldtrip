function [raw, event] = imotions2fieldtrip(filename, varargin)

% IMOTIONS2FIELDTRIP imports an iMotions *.txt file and represents it as a FieldTrip
% raw data structure.
%
% Use as
%   data = imotions2fieldtrip(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   fixtime       = 'no', 'squash' or 'interpolate' (default = 'no')
%   isnumeric     = cell-array with labels corresponding to numeric data (default = {})
%   isnotnumeric  = cell-array with labels not corresponding to numeric data (default = {})
%   isevent       = cell-array with labels corresponding to events (default = {})
%   isnotevent    = cell-array with labels not corresponding to events (default = {})
%
% The options 'isnumeric' and 'isnotnumeric' are mutually exclusive. Idem for
% 'isevent' and 'isnotevent'.
%
% When using fixtime='squash' identical timetamps are squashed together AND the whole
% timeaxis is interpolated to obtain a regularly sampled representation.
%
% See also FT_DATATYPE_RAW, FT_PREPROCESSING

% Copyright (C) 2017-2018, Robert Oostenveld
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

fixtime       = ft_getopt(varargin, 'fixtime', 'interpolate'); % squash, interpolate, none
isnotnumeric  = ft_getopt(varargin, 'isnotnumeric', {});
isnotevent    = ft_getopt(varargin, 'isnotevent', {});
isnumeric     = ft_getopt(varargin, 'isnumeric', {});
isevent       = ft_getopt(varargin, 'isevent', {});

% these options are mutually exclusive
if ~isempty(isnumeric) && ~isempty(isnotnumeric)
  error('you should specify either ''numeric'' or ''isnotnumeric''');
end
if ~isempty(isevent) && ~isempty(isnotevent)
  error('you should specify either ''isevent'' or ''isnotevent''');
end

% read the whole ASCII file into memory
% this will include a MATLAB table with the actual data
dat = read_imotions_txt(filename);

time       = dat.TimestampInSec;
label      = dat.table.Properties.VariableNames;
numericdat = zeros(0,numel(time));
numericsel = false(size(label));

if ~isempty(isnumeric)
  isnotnumeric = setdiff(label, isnumeric);
elseif ~isempty(isnotnumeric)
  isnumeric = setdiff(label, isnotnumeric);
end

if ~isempty(isevent)
  isnotevent = setdiff(label, isevent);
elseif ~isempty(isnotnumeric)
  isevent = setdiff(label, isnotevent);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for each field/column whether it is numerical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:numel(label)
  % skip if it is known to be not numeric
  if ismember(label{i}, isnotnumeric)
    continue
  end
  
  % don't convert if all empty
  str = dat.table.(label{i});
  if all(cellfun(@isempty, str))
    ft_info('column %15s does not contain numeric data', label{i});
    continue
  end
  
  % try converting the first element
  str = dat.table.(label{i})(1);
  val = str2double(str);
  if any(~cellfun(@isempty, str) & isnan(val))
    ft_info('column %15s does not contain numeric data', label{i});
    continue
  end
  
  % try converting the first 20 elements
  if numel(time)>20
    str = dat.table.(label{i})(1:20);
    val = str2double(str);
    if any(~cellfun(@isempty, str) & isnan(val))
      ft_info('column %15s does not contain numeric data', label{i});
      continue
    end
  end
  
  % try converting the whole column
  str = dat.table.(label{i});
  val = str2double(str);
  if all(cellfun(@isempty, str) | isnan(val))
    ft_info('column %15s does not contain numeric data', label{i});
    continue
  end
  
  % if it gets here, it means that the whole column is numerical
  numericsel(i) = true;
  ft_info('column %15s will be represented as channel', label{i});
  numericdat = cat(1, numericdat, val');
end
% remember the labels for the columns with numerical data
numericlabel = label(numericsel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct numerical channels for the columns that represent events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eventcode   = zeros(0,numel(time));
eventtype   = {};
eventvalue  = {};

% determine which channels are to be considered for events
eventsel = ~numericsel;
eventsel(ismember(label, isnotevent))   = false;
eventsel(strcmp(label, 'Timestamp'))    = false;
eventsel(strcmp(label, 'TimestampUTC')) = false;

for i=find(eventsel)
  str = dat.table.(label{i});
  if all(cellfun(@isempty, str))
    eventsel(i) = false;
    continue
  end
  
  % add one numerical channel per event type
  eventcode(end+1,:) = 0;
  eventtype{end+1}   = label{i};
  eventvalue{end+1}  = {};
  
  this = 1;
  code = 1; % this is the numerical code for the event values
  while this<=numel(str)
    
    next = find(~strcmp(str(this:end), str{this}), 1, 'first') + this - 1;
    if isempty(next)
      next = numel(str)+1;
    end
    
    % store the event as string and as numerical code
    eventvalue{end}{end+1} = str{this};
    eventcode(end,this:next-1) = code;
    
    this = next;
    code = code + 1;
  end
end

% give some feedback
for i=1:numel(eventtype)
  n = numel(eventvalue{i});
  if n>20
    % only give the summary
    ft_info('column %15s contains %d events, which are not shown in detail\n', eventtype{i}, n);
  else
    % give the full details
    ft_info('column %15s contains the following %d events\n', eventtype{i}, n);
    for j=1:numel(eventvalue{i})
      ft_info('%2d  %15s\n', j, eventvalue{i}{j});
    end
  end
end

% construct a raw data structure
raw.time{1}  = time(:)';
raw.trial{1} = cat(1, numericdat, eventcode);
raw.label    = cat(1, numericlabel(:), eventtype(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the next section deals with timestamps that are repeated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch fixtime
  case 'no'
    % keep it as it is
    
  case 'squash'
    % make a local copy for convenience
    time  = raw.time{1};
    trial = raw.trial{1};
    
    % the same timestamp can be on multiple lines in the file
    dt = diff(sort(time));
    
    if any(dt==0)
      ft_notice('squashing overlapping samples...\n');
      t = 1;
      while t<numel(time)
        sel = find(time==time(t));
        trial(:,t) = nanmean(trial(:,sel),2);
        time(sel(2:end)) = nan;
        t = t+numel(sel);
      end
      
      seltime = ~isnan(time);
      ft_notice('keeping %.0f %% of the original samples\n', 100*mean(seltime));
      time  = time (  seltime);
      trial = trial(:,seltime);
    end
    
    % put them back
    raw.time{1}  = time;
    raw.trial{1} = trial;
    
  case 'interpolate'
    y = raw.time{1};
    x = 1:numel(y);
    % use a GLM to estimate y = b0 * x + b1
    p = polyfit(x, y, 1);
    y = polyval(p, x);
    % replace the time by the estimated linear interpolant
    raw.time{1} = y;
    
  otherwise
    error('unsupported option')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate the data to ensure a regular time axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(fixtime, 'squash')
  % make a local copy for convenience
  time = raw.time{1};
  
  dt = diff(time);
  dt = median(dt);
  begtime = min(time);
  endtime = max(time);
  
  if any(diff(time)~=dt)
    ft_notice('resampling onto regularly spaced time axis\n');
    tmpcfg = [];
    tmpcfg.time = {begtime:dt:endtime};
    raw = ft_resampledata(tmpcfg, raw);
    [~, raw] = rollback_provenance([], raw);
  end
  
  % the channels with the event codes should remain integers
  sel = match_str(raw.label, eventtype);
  raw.trial{1}(sel,:) = floor(raw.trial{1}(sel,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct a structure with all events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

event = [];
nsample = numel(raw.time{1});
for i=1:numel(eventtype)
  eventcode = raw.trial{1}(i+numel(numericlabel),:);
  sel = [find(diff([0 eventcode])) nsample+1];
  for j=1:numel(sel)-1
    event(end+1).type     = eventtype{i};
    event(end  ).value    = eventvalue{i}{j};
    event(end  ).sample   = sel(j);
    event(end  ).duration = sel(j+1)-sel(j);
    event(end  ).offset   = 0;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wrap up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw.fsample = 1/median(diff(raw.time{1}));

% keep the details of the original tabular data
raw.hdr.orig = rmfield(dat, 'table');

% remove the channels with the integer representation of the events
raw.label    = raw.label(1:numel(numericlabel));
raw.trial{1} = raw.trial{1}(1:numel(numericlabel), :);
