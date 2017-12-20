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
  if ~isempty(str) && isnan(val)
    ft_info('column %s is not numeric', label{i});
    continue
  end
  
  % try converting the first 20 elements
  if numel(time)>10
    str = dat.data.(label{i})(1:20);
    val = str2double(str);
    if ~isempty(str) && any(isnan(val))
      ft_info('column %s is not numeric', label{i});
      continue
    end
  end
  
  % try converting the whole column
  str = dat.data.(label{i});
  val = str2double(str);
  if ~isempty(str) && any(isnan(val))
    ft_info('column %s is not numeric', label{i});
    continue
  end
  
  % if it gets here, it means that the whole column is numeric
  sellab(i) = true;
  ft_info('%s is numeric and will be represented as channel', label{i});
  numeric = cat(1, numeric, val');
end

% the same timestamp can be on multiple lines in the file
dt = diff(sort(time));

if any(dt==0)
  ft_notice('removing overlapping samples...\n');
  t = 1;
  while t<numel(time)
    sel = find(time==time(t));
    numeric(:,t) = nanmean(numeric(:,sel),2);
    time(sel(2:end)) = nan;
    t = t+numel(sel);
  end
  
  ft_notice('keeping %.0f of the original samples\n', mean(~isnan(time)));
  time    = time   (  ~isnan(time));
  numeric = numeric(:,~isnan(time));
end

% construct a raw data structure
raw.time = {time};
raw.trial = {numeric};
raw.label = label(sellab);

% interpolate the data
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
