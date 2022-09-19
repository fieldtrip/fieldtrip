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

% according to the eyelink documentation, the sample lines have different
% flavours, depending on how the data was collected (monocular or
% binocular), and possibly the conversion settings from edf2asc
% Monocular,                              <time> <xp>  <yp>  <ps>
% Monocular, with velocity                <time> <xp>  <yp>  <ps>  <xv>  <yv>
% Monocular, with resolution              <time> <xp>  <yp>  <ps>  <xr>  <yr>
% Monocular, with velocity and resolution <time> <xp>  <yp>  <ps>  <xv>  <yv>  <xr>  <yr>
% Binocular,                              <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr>
% Binocular, with velocity                <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr> <xvl> <yvl> <xvr> <yvr>
% Binocular, with resolution              <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr> <xr> <yr>
% Binocular, with velocity and resolution <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr> <xvl> <yvl> <xvr> <yvr><xr> <yr>
%
% Then, there may be additional columns, if recorded with cornea reflection
% mode:
%
% MONOCULAR Corneal Reflection (CR) Samples
% "..." if no warning for sample
% first character is "I" if sample was interpolated
% second character is "C" if CR missing
% third character is "R" if CR recovery in progress
%
% BINOCULAR Corneal Reflection (CR) Samples?
% "....." if no warning for sample
% first character is "I" if sample was interpolated
% second character is "C" if LEFT CR missing
% third character is "R" if LEFT CR recovery in progress
% fourth character is "C" if RIGHT CR missing
% fifth character is "R" if RIGHT CR recovery in progress
%
% Then, there may be additional columns, if data collection was done in
% remote mode:
%
% Data files recorded using the Remote Mode have extra columns to encode the
% target distance, position, and eye/target status information. The first three
% columns are:
% <target x>: X position of the target in camera coordinate (a value from 0 to 10000).
% Returns "MISSING_DATA" (-32768) if target is missing.
% <target y>: Y position of the target in camera coordinate (a value from 0 to 10000).
% Returns "MISSING_DATA" (-32768) if target is missing.
% <target distance>: Distance between the target and camera (in millimeters).
% Returns "MISSING_DATA" (-32768) if target is missing.
% The next thirteen fields represent warning messages for that sample relating to
% the target and eye image processing."............." if no warning for target and eye image
% first character is "M" if target is missing
% second character is "A" if extreme target angle occurs
% third character is "N" if target is near eye so that the target window and eye window overlap
% fourth character is "C" if target is too close
% fifth character is "F" if target is too far
% sixth character is "T" if target is near top edge of the camera image
% seventh character is "B" if target is near bottom edge of the camera image
% eighth character is "L" if target is near left edge of the camera image
% ninth character is "R" if target is near right edge of the camera image
% tenth character is "T" if eye is near top edge of the camera image
% eleventh character is "B" if eye is near bottom edge of the camera image
% twelfth character is "L" if eye is near left edge of the camera image
% thirteenth character is "R" if eye is near right edge of the camera image
% For a binocular recording, there will be seventeen target/eye status columns,
% with the last eight columns reporting the warning messages for the left and
% right eyes separately.

% So, anecdotally, parsing the file is a bit tedious, and not very robust.
% Reading and converting per line is super slow, so here I chose to pipe
% the datalines through a checker that removes all '...' '.....'
% '.............' and '.................' (or containing warning) stuff
% so that we end up with more or less well behaved data.


% check if the formatting for all datlines is similar
seltab = char(datline)==sprintf('\t');
ntab   = sum(seltab,2);

% check whether all lines have the same number of columns
assert(all(ntab==ntab(1)));

% identify the chunks of consecutive '...' and comments
str = char(datline)';
str(end+1, :) = sprintf('\t');
siz = size(str);
str = uint8(str(:)');

boolval = str==46 | (str>=65&str<=90); % dots or capital letters
begsmp  = find(boolval==1 & [0 boolval(1:end-1)]==0);
endsmp  = find(boolval==1 & [boolval(2:end) 1]==0);
sel     = ~ismember(endsmp-begsmp+1, [3 5 13 17]);
boolval(begsmp(sel)) = false;
str(boolval) = 32;

str = reshape(char(str), siz); % now all the warnings etc should have been removed
datline = cellstr(str(1:end-1,:)');

% check again if the formatting for all datlines is similar
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

C   = textscan(datline, fmt, nline, 'Delimiter', sprintf('\t'));
dat = zeros(nline, numel(C)) + nan;
for i=1:numel(C)
  notok = strcmp(C{i}, '.')|strcmp(C{i}, '  ')|strcmp(C{i}, '     ')|cellfun('isempty', C{i})|startsWith(C{i}, '. ');
  ok    = ~notok;
  tmp = reshape([char(C{i}(ok))';repmat(' ',1,sum(ok))],[],1)';

  % convert to floats, much faster than str2double
  if sum(ok)
    tmp = textscan(tmp, '%f');
    dat(ok, i) = tmp{1};
  end
end
asc.dat = dat';

% remove any all-nan rows
asc.dat(sum(~isfinite(asc.dat),2)==size(asc.dat,2), :) = [];

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
