function dat = read_imotions_txt(filename)

% READ_IMOTIONS_TXT reads *.txt files that are exported from the iMotions software.
%
% Use as
%   dat = read_imotions_txt(filename
%
% See also TEXTSCAN

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

fprintf('reading header from %s\n', filename);
% read the whole file at once
fid = fopen(filename, 'rt');
aline = fread(fid, inf, 'char=>char');          % returns a single long string
fclose(fid);

aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
aline = tokenize(aline, uint8(newline));        % split on newline

dat = [];

% the first section contains multiple lines with general header information
headerline = 1;
while startsWith(aline{headerline}, '#')
  % for example "#Version : 6.3.6059.1"
  tok = tokenize(aline{headerline}(2:end), ':');
  key = matlab.lang.makeValidName(strtrim(tok{1}));
  val = strtrim(tok{2});
  dat.(key) = val;
  headerline = headerline + 1;
end

% there is one line with column headings
headerline = 1;
while ~startsWith(aline{headerline}, 'StudyName')
  headerline = headerline+1;
end
dat.Labels = tokenize(aline{headerline}, uint8(9)); % split on tab, which is ascii code 9
dat.Labels = matlab.lang.makeValidName(dat.Labels);

fprintf('parsing tabular data\n');
% now that the header information is known, the actual content can be scanned
fid = fopen(filename, 'rt');
c = textscan(fid, repmat('%s', [1 numel(dat.Labels)]), 'Delimiter', '\t', 'Headerlines', headerline);
fclose(fid);

% add each column as a variable to the data structure
for i=1:numel(dat.Labels)
  dat.table.(dat.Labels{i}) = c{i};
end
% convert to a table object (requires MATLAB R2013b and up)
dat.table = struct2table(dat.table);

fprintf('converting timestamps to seconds\n');
% convert timestamp to seconds
% the first call is not accurate enough to represent the miliseconds
% dat.TimestampInSec = datenum(dat.table.Timestamp, 'yyyymmdd_HHMMSSFFF');
time = datevec(dat.table.Timestamp, 'yyyymmdd_HHMMSSFFF');
time = time(:,4)*60*60 + time(:,5)*60 + time(:,6);
dat.TimestampInSec = time - time(1);




