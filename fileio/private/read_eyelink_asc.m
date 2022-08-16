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

% read the whole file at once
fid = fopen_or_error(filename, 'rt');
aline = fread(fid, inf, 'char=>char');          % returns a single long string
fclose(fid);

aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
aline = tokenize(aline, uint8(newline));        % split on newline

sel = startsWith(aline, '0') | ...
      startsWith(aline, '1') | ...
      startsWith(aline, '2') | ...
      startsWith(aline, '3') | ...
      startsWith(aline, '4') | ...
      startsWith(aline, '5') | ...
      startsWith(aline, '6') | ...
      startsWith(aline, '7') | ...
      startsWith(aline, '8') | ...
      startsWith(aline, '9');

datline = aline(sel);
aline   = aline(~sel);

% check if the formatting for all datlines is similar
nline  = numel(datline);
seltab = char(datline)==sprintf('\t');
ntab   = sum(seltab,2);

% check whether all lines have the same number of columns
assert(all(ntab==ntab(1)));

fmt = repmat('%s ', [1 ntab(1)+1]);
fmt = fmt(1:end-1);

% convert datline in a long string again
datline = char(datline)';
datline(end+1, :) = sprintf('\t');
datline = datline(:)';

C   = textscan(datline, fmt, nline);
dat = zeros(nline, numel(C)) + nan;
for i=1:numel(C)
  notok = strcmp(C{i}, '.')|strcmp(C{i}, '...');
  tmp = reshape([char(C{i}(~notok))';repmat(' ',1,sum(~notok))],[],1)';

  % convert to floats, much faster than str2double
  if ~isempty(tmp)
    tmp = textscan(tmp, '%f');
    dat(~notok, i) = tmp{1};
  end
end
asc.dat = dat';

selsfix    = startsWith(aline, 'SFIX');
selefix    = startsWith(aline, 'EFIX');
selssacc   = startsWith(aline, 'SSACC');
selesacc   = startsWith(aline, 'ESACC');
selsblink  = startsWith(aline, 'SBLINK');
seleblink  = startsWith(aline, 'EBLINK');
selmsg     = startsWith(aline, 'MSG');
selhdr     = startsWith(aline, '**');
selinput   = startsWith(aline, 'INPUT');

asc.sfix   = tocell(aline(selsfix)); % convert cell-vector of lines into a cell array
asc.efix   = tocell(aline(selefix));
asc.ssacc  = tocell(aline(selssacc));
asc.esacc  = tocell(aline(selesacc));
asc.sblink = tocell(aline(selsblink));
asc.eblink = tocell(aline(seleblink));
asc.msg    = tocell(aline(selmsg), 1);
asc.header = char(aline(selhdr));
asc.input  = tocell(aline(selinput));

% convert the cell-arrays into tables
if ~isempty(asc.sfix),   asc.sfix   = totable(asc.sfix(:, 2:end),   {'eye' 'stime'}); end
if ~isempty(asc.sblink), asc.sblink = totable(asc.sblink(:, 2:end), {'eye' 'stime'}); end
if ~isempty(asc.ssacc),  asc.ssacc  = totable(asc.ssacc(:, 2:end),  {'eye' 'stime'}); end
if ~isempty(asc.input),  asc.input  = totable(asc.input(:, 2:end),  {'timestamp', 'value'}, [1 2]); end

% these lines may have been formatted in different ways
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

% remove any all-nan rows
asc.dat(sum(~isfinite(asc.dat),2)==size(asc.dat,2), :) = [];

function lineout = convertline(linein, msgflag)

% split the line by horizontal tabs, and the first output element
% thereof by spaces. 

%tmp     = strrep(linein, sprintf('\t'), ' ');
lineout = strsplit(linein, '\t');
lineout = [strsplit(lineout{1}, ' ') lineout(2:end)];
if msgflag
  tmp = strsplit(lineout{2}, ' ');
  lineout = [lineout(1) tmp(1) {sprintf('%s ', tmp{2:end})}];
end

function tableout = totable(cellin, fnames, convert)

if nargin<3
  convert = 2:numel(fnames);
end

% convert into table, and do a str2double conversion for the numeric data,
% assuming only the first variable to be kept as a string

tableout = cell2table(cellin, 'VariableNames', fnames);
for k = convert
  tableout.(fnames{k}) = str2double(tableout.(fnames{k}));
end

function cellout = tocell(textin, msgflag)

if nargin<2
  msgflag = false;
end

cellout = {};
for k = 1:numel(textin)
  tline   = convertline(textin{k}, msgflag); 
  if k==1
    cellout = cell(numel(textin), numel(tline));
  end
  cellout(k,:) = tline;
end
