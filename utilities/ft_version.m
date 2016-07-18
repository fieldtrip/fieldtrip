function [ftver, ftpath] = ft_version

% FT_VERSION returns the version and installation directory of FieldTrip
%
% FieldTrip is not released with version numbers as "2.0", "2.1", etc. Instead, we
% have a Subversion (SVN) development version and a daily FTP release version.
%
% The SVN development version is labeled with the revision number like "rXXXXX",
% where XXXX is the revision number.
%
% The daily FTP release version is packaged as a zip file and its version is
% indicated with "YYMMDD" (year, month, day).
%
% Use as
%   ft_version
% to display the latest revision number on screen, or
%   [ftver, ftpath] = ft_version
% to get the version and the installation root directory.
%
% See also VERSION, VER

% Copyright (C) 2012, Eelke Spaak
%
% This file is part of FieldTrip, see http://www.ru.nl/donders/fieldtrip
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

persistent issvn
persistent isgit

ftpath = fileparts(mfilename('fullpath'));
ftpath = ftpath(1:end-10); % strip away '/utilities' where this function is located

if isempty(issvn)
  % are we dealing with an SVN working copy of fieldtrip?
  issvn = isdir(fullfile(ftpath, '.svn'));
end

if isempty(isgit)
  % are we dealing with an GIT working copy of fieldtrip?
  isgit = isdir(fullfile(ftpath, '.git'));
end

if ispc
  % this requires a file extension
  ext = '.exe';
else
  ext = '';
end

if issvn
  % use svn system call to determine latest revision
  olddir = pwd();
  cd(ftpath);
  [status, output] = system(sprintf('svn%s info', ext));
  cd(olddir);
  if status > 0
    if ~ispc
      % the command line tools will probably not be available on windows
      warning('you seem to have an SVN development copy of FieldTrip, yet ''svn info'' does not work as expected');
    end
    ftver = 'unknown';
  else
    rev = regexp(output, 'Revision: (.*)', 'tokens', 'dotexceptnewline');
    rev = rev{1}{1};
    ftver = ['r' rev];
  end
  
elseif isgit
  % use git system call to determine latest revision
  olddir = pwd();
  cd(ftpath);
  [status, output] = system(sprintf('git%s rev-parse --short HEAD', ext));
  cd(olddir);
  if status > 0
    if ~ispc
      % the command line tools will probably not be available on windows
      warning('you seem to have an GIT development copy of FieldTrip, yet ''git rev-parse'' does not work as expected');
    end
    ftver = 'unknown';
  else
    ftver = strtrim(output); % remove trailing newline character
  end
  
elseif isequal(regexp(ftpath, ['.*' filesep 'fieldtrip-fieldtrip-[[0-9][a-z]]{7}']), 1)
  % this corresponds with being downloaded from the Mathworks file exchange link to github
  % which results in a ftpath like /Users/robert/matlab/fieldtrip-fieldtrip-851478d
  ftver = ftpath(end-6:end);
  
elseif isequal(regexp(ftpath, ['.*' filesep 'fieldtrip-20[0-9]{6}']), 1)
  % this corresponds with the daily version from the ftp server
  % which results in a ftpath like /Users/robert/matlab/fieldtrip-20160317
  ftver = ftpath(end-7:end);
  
else
  % get it from the Contents.m file in the FieldTrip release
  a = ver(ftpath);
  ftver = a.Version;
  
end % if issvn, isgit or otherwise

if nargout==0
  fprintf('\nThis is FieldTrip, version %s.\n\n', ftver);
  clear ftver
end
