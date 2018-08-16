function event = presentation_log(filename)

% PRESENTATION_LOG reads a NBS Presentation scenario log file and
% represents it as a FieldTrip event structure.
%
% See also FT_READ_EVENT

% Copyright (C) 2018 Robert Oostenveld
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

fid = fopen(filename, 'rt');

line = {};
while ~feof(fid)
  line = cat(1, line, fgetl(fid));
end
fclose(fid);

% remove empty lines
sel = cellfun(@isempty, line);
line = line(~sel);

header1 = find(startsWith(line, 'Subject'));
header2 = find(startsWith(line, 'Event'));

part1 = line((header1+1):(header2-1));
part2 = line((header2+1):(end));

VariableNames = tokenize(line{header1}, sprintf('\t'));
for i=1:length(VariableNames)
  VariableNames{i} = fixname(VariableNames{i});
end

Values = cell(length(part1), length(VariableNames));
for i=1:length(part1)
  val = tokenize(part1{i}, sprintf('\t'));
  val((end+1):length(VariableNames)) = {nan};
  Values(i,:) = val;
end

table1 = cell2table(Values);
table1.Properties.VariableNames = uniquelabels(VariableNames);

% where possible convert columns from string into numeric values
for i=1:numel(VariableNames)
  numeric = cellfun(@str2double, table1.(VariableNames{i}));
  if ~any(isnan(numeric))
    table1.(VariableNames{i}) = numeric;
  end
end

numevent  = size(table1,1);
type      = table1.event_type;
value     = table1.code;
sample    = num2cell(nan(numevent, 1));
timestamp = num2cell(table1.time);
offset    = cell(numevent, 1);
duration  = cell(numevent, 1);

event = struct('type', type, 'value', value, 'sample', sample, 'timestamp', timestamp, 'offset', offset, 'duration', duration);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = uniquelabels(c)
for i=1:numel(c)-1
  sel = find(strcmp(c((i+1):end), c{i}));
  if ~isempty(sel)
    for j=1:numel(sel)
      c{sel(j)+i} = sprintf('%s_%d', c{sel(j)+i}, j+1);
    end
  end
end
