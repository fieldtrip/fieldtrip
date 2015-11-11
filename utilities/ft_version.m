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

% show the latest revision present in this copy of fieldtrip

if issvn
  % use svn system call to determine latest revision
  olddir = pwd();
  cd(ftpath);
  [status, output] = system('svn info');
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
  tmpfile = tempname;
  
  olddir = pwd();
  cd(ftpath);
  [status, output] = system(sprintf('git show > %s', tmpfile));
  cd(olddir);
  if status > 0
    % FIXME the command line tools will probably not be available on windows
    error('you seem to have an GIT development copy of FieldTrip, yet ''git show'' does not work as expected');
  end
  
  fp = fopen(tmpfile);
  if fp>0
    line = fgetl(fp); % the first line contains the commit number
    fclose(fp);
    rev = regexp(line, ' ', 'split');
    rev = rev{2};
    
    % this is a string like 4d3c309129f12146885120c2853a11362e048ea7
    ftver = rev;
  else
    ftver = 'unknown';
  end
  
else
  % get it from the Contents.m file in the FieldTrip release
  a = ver(ftpath);
  ftver = a.Version;
  
end % if issvn, isgit or otherwise

if nargout==0
  fprintf('\nThis is FieldTrip, version %s.\n\n', ftver);
  clear ftver
end


