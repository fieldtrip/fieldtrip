function asc = read_eyelink_asc(filename)

% READ_EYELINK_ASC reads the header information, input triggers, messages
% and all data points from an Eyelink *.asc file. The output events are
% represented as matlab tables (after Aug 2022)
%
% Use as
%   asc = read_eyelink_asc(filename)

% Copyright (C) 2010-2015, Robert Oostenveld
% Copyright (C) 2022, Jan-Mathijs Schoffelen
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

asc.header  = {};
asc.msg     = {};
asc.input   = [];
asc.sfix    = {};
asc.efix    = {};
asc.ssacc   = {};
asc.esacc   = {};
asc.dat     = [];
asc.sblink   = {}; % blink parsing added
asc.eblink   = {};
current   = 0;

% read the whole file at once
fid = fopen_or_error(filename, 'rt');
aline = fread(fid, inf, 'char=>char');          % returns a single long string
fclose(fid);

aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
aline = tokenize(aline, uint8(newline));        % split on newline

for i=1:numel(aline)
  % this is a bit inefficient, due to the massive cat and sscanf operations, but I
  % keep it like this for now

  tline = aline{i};
  
  if numel(tline) && any(tline(1)=='0':'9')
    % if regexp(tline, '^[0-9]')
    tline   = strrep(tline, ' . ', ' NaN '); % replace missing values
    tmp     = sscanf(tline, '%f');
    nchan   = numel(tmp);
    current = current + 1;
    
    if size(asc.dat, 1)<nchan
      % increase the allocated number of channels
      asc.dat(nchan,:) = 0;
    end
    
    if size(asc.dat, 2)<current
      % increase the allocated number of samples
      asc.dat(:,end+10000) = 0;
    end
    
    % add the current sample to the data matrix
    asc.dat(1:nchan, current) = tmp;
    
    
  elseif regexp(tline, '^INPUT')
    [val, num] = sscanf(tline, 'INPUT %d %d');
    this.timestamp = val(1);
    this.value     = val(2);
    if isempty(asc.input)
      asc.input = this;
    else
      asc.input = cat(1, asc.input, this);
    end
    
  elseif regexp(tline, '\*\*.*')
    asc.header = cat(1, asc.header, {tline});
    
  elseif regexp(tline, '^MSG')
    asc.msg = cat(1, asc.msg, {tline});
    
  elseif regexp(tline, '^SFIX')
    tline    = convertline(tline);
    asc.sfix = cat(1, asc.sfix, tline);
    
  elseif regexp(tline, '^EFIX')
    tline    = convertline(tline);
    asc.efix = cat(1, asc.efix, tline);
    
  elseif regexp(tline, '^SSACC')
    tline     = convertline(tline);
    asc.ssacc = cat(1, asc.ssacc, tline);
    
  elseif regexp(tline, '^ESACC')
    tline     = convertline(tline);
    asc.esacc = cat(1, asc.esacc, tline);
    
  elseif regexp(tline, '^SBLINK')
    tline      = convertline(tline);
    asc.sblink = cat(1, asc.sblink, tline);
    
  elseif regexp(tline, '^EBLINK')
    tline      = convertline(tline);
    asc.eblink = cat(1, asc.eblink, tline);
    
  else
    % all other lines are not parsed
  end
  
end

if ~isempty(asc.input)
  asc.input = struct2table(asc.input);
end

% convert the cell-arrays into tables
if ~isempty(asc.sfix),   asc.sfix   = totable(asc.sfix(:, 2:end),   {'eye' 'stime'}); end
if ~isempty(asc.sblink), asc.sblink = totable(asc.sblink(:, 2:end), {'eye' 'stime'}); end
if ~isempty(asc.ssacc),  asc.ssacc  = totable(asc.ssacc(:, 2:end),  {'eye' 'stime'}); end

if ~isempty(asc.efix)
  try
    asc.efix = totable(asc.efix(:,2:end), {'eye', 'stime', 'etime', 'dur', 'axp', 'ayp', 'aps'});
  catch
    asc.efix = totable(asc.efix(:,2:end), {'eye', 'stime', 'etime', 'dur', 'axp', 'ayp', 'aps', 'xr', 'yr'});
  end
end
if ~isempty(asc.eblink)
  asc.eblink = totable(asc.eblink(:,2:end), {'eye', 'stime', 'etime', 'dur'});
end
if ~isempty(asc.esacc)
  try
    asc.esacc = totable(asc.esacc(:,2:end), {'eye', 'stime', 'etime', 'dur', 'sxp', 'syp', 'exp', 'eyp', 'ampl', 'pv'});
  catch
    asc.esacc = totable(asc.esacc(:,2:end), {'eye', 'stime', 'etime', 'dur', 'sxp', 'syp', 'exp', 'eyp', 'ampl', 'pv', 'xr', 'yr'});
  end
end

% remove the samples that were not filled with real data
asc.dat = asc.dat(:,1:current);

function lineout = convertline(linein)

% split the line by horizontal tabs, and the first output element
% thereof by spaces. 

tmp     = strrep(linein, sprintf('\t'), ' ');
lineout = strsplit(tmp, ' ');
lineout = lineout(:)';

function tableout = totable(cellin, fnames)

% convert into table, and do a str2double conversion for the numeric data,
% assuming only the first variable to be kept as a string

tableout = cell2table(cellin, 'VariableNames', fnames);
for k = 2:numel(fnames)
  tableout.(fnames{k}) = str2double(tableout.(fnames{k}));
end
