function [ftver, ftpath] = ft_version(command)

% FT_VERSION returns the version and installation directory of FieldTrip
%
% FieldTrip is not released with version numbers as "2.0", "2.1", etc. Instead, we
% share our development version on http://github.com/fieldtrip. You can use git or
% subversion (svn) to make a local version of the repository. Furthermore, we release
% daily version as zip-file on our FTP server.
%
% If you access the development version using git, it is labeled with the hash of the
% latest commit like "128c693". You can access the specific version "XXXXXX" at
% https://github.com/fieldtrip/fieldtrip/commit/XXXXXX.
%
% If you access the development version using svn, it is labeled with the revision
% number like "rXXXXX", where XXXX is the revision number.
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
% When using git for version control, you can also get additional information with
%   ft_version revision
%   ft_version branch
%   ft_version clean
%
% See also FT_PLATFORM_SUPPORTS, VERSION, VER, VERLESSTHAN

% Copyright (C) 2012-2016, Eelke Spaak
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

if nargin<1
  % this is only supported for git
  command='revision';
end

ftpath = fileparts(mfilename('fullpath'));
ftpath = ftpath(1:end-10); % strip away '/utilities' where this function is located

if isempty(issvn)
  % are we dealing with an SVN working copy of FieldTrip?
  issvn = isfolder(fullfile(ftpath, '.svn'));
end

if isempty(isgit)
  % are we dealing with an GIT working copy of FieldTrip?
  isgit = exist(fullfile(ftpath, '.git'), 'file');
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
      ft_warning('you seem to have an SVN development copy of FieldTrip, yet ''svn info'' does not work as expected');
    end
    ftver = 'unknown';
  else
    rev = regexp(output, 'Revision: (.*)', 'tokens', 'dotexceptnewline');
    rev = rev{1}{1};
    ftver = ['r' rev];
  end
  
elseif isgit
  % test whether the git executable is available
  [status, output] = system(sprintf('git%s --version', ext));
  if status>0
    if ~ispc
      % the command line tools will probably not be available on windows
      ft_warning('you seem to have an GIT development copy of FieldTrip, yet ''git'' does not work as expected');
    end
    ftver = 'unknown';
    
  else
    % use git system call to determine latest revision
    olddir = pwd();
    cd(ftpath);
    switch command
      case 'branch'
        [status, output] = system(sprintf('git%s rev-parse --abbrev-ref HEAD', ext));
        ftver = strtrim(output); % remove trailing newline character
      case 'revision'
        [status, output] = system(sprintf('git%s rev-parse --short HEAD', ext));
        ftver = strtrim(output); % remove trailing newline character
      case 'clean'
        [status, output] = system(sprintf('git%s diff --quiet --exit-code', ext));
        if status
          ftver = 'no';
        else
          ftver = 'yes';
        end
      otherwise
        ft_error('unsupported command "%s"');
    end
    cd(olddir);
    
  end % if git available
  
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
  fprintf('\nThis is FieldTrip, %s %s.\n\n', command, ftver);
  clear ftver ftpath
end
