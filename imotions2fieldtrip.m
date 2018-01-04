function [raw, event] = imotions2fieldtrip(filename, varargin)

% IMOTIONS2FIELDTRIP imports an iMotions *.txt file and represents it as a FieldTrip
% raw data structure.
%
% Use as
%   data = imotions2fieldtrip(filename)
%
% See also FT_DATATYPE_RAW, FT_PREPROCESSING

% Copyright (C) 2017, Robert Oostenveld
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

% read the whole ASCII file into memory
% this will include a MATLAB table with the actual data
dat = read_imotions_txt(filename);

time    = dat.TimestampInSec;
numeric = zeros(0,numel(time));
label   = dat.data.Properties.VariableNames;
sellab  = false(size(label));

% check for each field/column whether it is numeric
for i=1:numel(label)
  % try converting the first element
  str = dat.data.(label{i})(1);
  val = str2double(str);
  if any(~cellfun(@isempty, str) & isnan(val))
    ft_info('column %15s does not contain numeric data', label{i});
    continue
  end
  
  % try converting the first 20 elements
  if numel(time)>10
    str = dat.data.(label{i})(1:20);
    val = str2double(str);
    if any(~cellfun(@isempty, str) & isnan(val))
      ft_info('column %15s does not contain numeric data', label{i});
      continue
    end
  end
  
  % try converting the whole column
  str = dat.data.(label{i});
  val = str2double(str);
  if all(cellfun(@isempty, str) | isnan(val))
    ft_info('column %15s does not contain numeric data', label{i});
    continue
  end
  
  % if it gets here, it means that the whole column is numeric
  sellab(i) = true;
  ft_info('column %15s will be represented as channel', label{i});
  numeric = cat(1, numeric, val');
end
% remember the labels for the columns with numeric data
numericlabel = label(sellab);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct numeric channels for the columns that represent events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eventcode   = zeros(0,numel(time));
eventtype   = {};
eventvalue  = {};

% these are not to be considered for events
sellab(strcmp(label,'Timestamp'))    = true;
sellab(strcmp(label,'UTCTimestamp')) = true;

for i=find(~sellab)
  str = dat.data.(label{i});
  if all(cellfun(@isempty, str))
    continue
  end
  
  % add one numeric channel per event type
  eventcode(end+1,:) = 0;
  eventtype{end+1}   = label{i};
  eventvalue{end+1}  = {};
  
  this = 1;
  code = 1; % this is the numeric code for the event values
  while this<numel(str)
    
    next = find(~strcmp(str(this:end), str{this}), 1, 'first') + this;
    if isempty(next)
      next = numel(str)+1;
    end
    
    % store the event as string and as numeric code
    eventvalue{end}{end+1} = str{this};
    eventcode(end,this:next-1) = code;
    
    this = next;
    code = code + 1;
  end
end

% construct a raw data structure
raw.time{1}  = time;
raw.trial{1} = cat(1, numeric, eventcode);
raw.label    = cat(1, numericlabel(:), eventtype(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the next section deals with timestamps that are repeated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a local copy for convenience
time  = raw.time{1};
trial = raw.trial{1};

% the same timestamp can be on multiple lines in the file
dt = diff(sort(time));

if any(dt==0)
  ft_notice('removing overlapping samples...\n');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate the data to ensure a regular time axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct a structure with all events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

event = [];
nsample = numel(time)+1;
for i=1:numel(eventtype)
  eventcode = raw.trial{1}(i+numel(numericlabel),:);
  sel = [find(diff([0 eventcode])) nsample+1];
  for j=1:numel(sel)-1
    event(end+1).type     = eventtype{i};
    event(end  ).value    = eventvalue{i}{j};
    event(end  ).sample   = sel(j);
    event(end  ).duration = sel(j+1)-sel(j)+1;
    event(end  ).offset   = 0;
  end
end

% the channels with the integer events should be removed, but for now I am
% keeping them to help with debugging

