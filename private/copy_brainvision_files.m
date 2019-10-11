function copy_brainvision_files(oldheaderfile, newheaderfile, deleteflag)

% COPY_BRAINVISION_FILES copies a BrainVision EEG dataset, which consists of a vhdr
% header file, vmrk marker file and a data file with the extension dat, eeg or seg.
%
% Use as
%   copy_brainvision_files(oldname, newname, deleteflag)
%
% Both the old and the new filename should be strings corresponding to the header
% file, i.e. including the vhdr extension.
%
% The third "deleteflag" argument is optional, it should be a boolean
% that specifies whether the original files should be deleted after
% copying or not (default = false).
%
% An earlier version of this function can be found on
%   - https://gist.github.com/robertoostenveld/e31637a777c514bf1e86272e1092316e
%   - https://gist.github.com/CPernet/e037df46e064ca83a49fb4c595d4566a
%
% See also COPY_CTF_FILES

% Copyright (C) 2018-2019, Robert Oostenveld
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

%% deal with inputs

if nargin<3
  deleteflag = false;
end

if ~endsWith(oldheaderfile, 'vhdr') && ~endsWith(oldheaderfile, 'VHDR')
  ft_error('you should specify the old header file');
end

if ~endsWith(newheaderfile, 'vhdr') && ~endsWith(newheaderfile, 'VHDR')
  ft_error('you should specify the new header file');
end

% determine whether the file extensions should be in lower or upper case
if ~isempty(regexp(newheaderfile, 'VHDR$', 'once'))
  switchcase = @upper;
else
  switchcase = @lower;
end

% determine the file name and path
[po, fo, xo] = fileparts(oldheaderfile);
[pn, fn, xn] = fileparts(newheaderfile);

%% copy the header file

if ~exist(oldheaderfile, 'file')
  ft_error('the file "%s" does not exists', oldheaderfile);
end

if exist(newheaderfile, 'file')
  ft_error('the file "%s" already exists', newheaderfile);
end

fid1 = fopen(oldheaderfile, 'r'); % read old
fid2 = fopen(newheaderfile, 'w'); % write new

while ~feof(fid1)
  line = fgetl(fid1);
  if ~isempty(regexp(line, '^MarkerFile', 'once'))
    oldmarkerfile = getvalue(line);
    [xx, xx, xo] = fileparts(oldmarkerfile); % ignore the name, just get the extension
    oldmarkerfile = [fo switchcase(xo)];
    newmarkerfile = [fn switchcase(xo)];
    line = sprintf('MarkerFile=%s', newmarkerfile);
  elseif ~isempty(regexp(line, '^DataFile', 'once'))
    olddatafile = getvalue(line);
    [xx, xx, xo] = fileparts(olddatafile); % ignore the name, just get the extension
    olddatafile = [fo switchcase(xo)];
    newdatafile = [fn switchcase(xo)];
    line = sprintf('DataFile=%s', newdatafile);
  end
  fprintf(fid2, '%s\r\n', line);
end
fclose(fid1);
fclose(fid2);

%% copy the marker file

% use the full path
oldmarkerfile = fullfile(po, oldmarkerfile);
newmarkerfile = fullfile(pn, newmarkerfile);

if ~exist(oldmarkerfile, 'file')
  ft_error('the file "%s" does not exist', oldmarkerfile);
end

if exist(newmarkerfile, 'file')
  ft_error('the file "%s" already exists', newmarkerfile);
end

fid1 = fopen(oldmarkerfile, 'r');
fid2 = fopen(newmarkerfile, 'w');

while ~feof(fid1)
  line = fgetl(fid1);
  if ~isempty(regexp(line, '^DataFile', 'once'))
    olddatafile = getvalue(line);
    [xx, xx, xo] = fileparts(olddatafile); % ignore the name, just get the extension
    olddatafile = [fo switchcase(xo)];
    newdatafile = [fn switchcase(xo)];
    line = sprintf('DataFile=%s', newdatafile);
  end
  fprintf(fid2, '%s\r\n', line);
end
fclose(fid1);
fclose(fid2);

%% copy the data file

% use the full path
olddatafile = fullfile(po, olddatafile);
newdatafile = fullfile(pn, newdatafile);

if ~exist(olddatafile, 'file')
  error('the file "%s" does not exists', olddatafile);
end

if exist(newdatafile, 'file')
  ft_error('the file "%s" already exists', newdatafile);
end

status = copyfile(olddatafile, newdatafile);

if ~status
  error('failed to copy data from "%s" to "%s"', olddatafile, newdatafile);
end

%% delete the old files
if deleteflag
  try
    delete(oldheaderfile);
  catch
    ft_warning('cannot delete "%s"', oldheaderfile);
  end
  try
    delete(oldmarkerfile);
  catch
    ft_warning('cannot delete "%s"', oldmarkerfile);
  end
  try
    delete(olddatafile);
  catch
    ft_warning('cannot delete "%s"', olddatafile);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION retuen the part after the '='
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = getvalue(str)
tok = tokenize(str, '=');
val = tok{2};
