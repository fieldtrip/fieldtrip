function varargout = ft_version(command)

% FT_VERSION returns the version of FieldTrip and the path where it is installed
%
% FieldTrip is not released with version numbers as "2.0", "2.1", etc. Instead, we
% share our development version on http://github.com/fieldtrip/fieldtrip. You can use
% git to make a local clone of the development version. Furthermore, we make
% more-or-less daily releases of the code available on
% https://github.com/fieldtrip/fieldtrip/releases and as zip file on our FTP server.
%
% If you use git with the development version, the version is labeled with the hash
% of the latest commit like "128c693". You can access the specific version "XXXXXX"
% at https://github.com/fieldtrip/fieldtrip/commit/XXXXXX.
%
% If you download the daily released version from our FTP server, the version is part
% of the file name "fieldtrip-YYYYMMDD.zip", where YYY, MM and DD correspond to year,
% month and day.
%
% Use as
%   ft_version
% to display the latest revision number on screen, or
%   [ftver, ftpath] = ft_version
% to get the version and the installation root directory.
%
% When using git with the development version, you can also get additional information with
%   ft_version revision
%   ft_version branch
%   ft_version clean
%
% See also FT_PLATFORM_SUPPORTS, VERSION, VER, VERLESSTHAN

% Copyright (C) 2012-2019, Eelke Spaak, Robert Oostenveld
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
persistent ftver
persistent ftpath
persistent is_pc

if isempty(is_pc)
  is_pc = ispc;
end

if isempty(ftpath)
  ftpath = fileparts(mfilename('fullpath'));
  ftpath = ftpath(1:end-10); % strip away '/utilities' where this function is located
end

% set the defaults
clean = 'no';
branch = 'unknown';

if is_pc
  % this requires a file extension
  ext = '.exe';
else
  ext = '';
end

if isempty(issvn)
  % are we dealing with an SVN working copy of FieldTrip?
  issvn = isfolder(fullfile(ftpath, '.svn'));
end

if isempty(isgit)
  % are we dealing with an GIT working copy of FieldTrip?
  isgit = exist(fullfile(ftpath, '.git'), 'file');
  % is the git command line executable available?
  if isgit
    [status, output] = system(sprintf('git%s --version', ext));
    if status>0
      if ~is_pc
        % the git command line executable will probably not be available on windows
        ft_warning('you seem to have an GIT development copy of FieldTrip, yet ''git'' does not work as expected');
      end
      isgit = false;
    end
  end
end

if nargin<1
  % this is only supported for git
  command = 'revision';
end

if ~isempty(ftver) && ~isempty(ftpath) && nargin<1
  % use the previously determined values
  
elseif issvn
  % use svn system call to determine latest revision
  [status, output] = system(sprintf('cd %s && svn%s info', ftpath, ext));
  if status > 0
    if ~is_pc
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
  % use git system call to determine latest revision
  switch command
    case 'revision'
      [status, output] = system(sprintf('cd %s && git%s rev-parse --short HEAD', ftpath, ext));
      ftver = strtrim(output); % remove trailing newline character
    case 'branch'
      [status, output] = system(sprintf('cd %s && git%s rev-parse --abbrev-ref HEAD', ftpath, ext));
      branch = strtrim(output); % remove trailing newline character
    case 'clean'
      [status, output] = system(sprintf('cd %s && git%s diff --quiet --exit-code', ftpath, ext));
      if status
        clean = 'no';
      else
        clean = 'yes';
      end
    otherwise
      ft_error('unsupported command "%s"', command);
  end % switch command
  
elseif isequal(regexp(ftpath, ['.*\' filesep '[fF]ieldtrip-fieldtrip-[[0-9][a-z]]{7}']), 1)
  % this corresponds with being downloaded from the Mathworks file exchange link to github
  % which results in a ftpath like /Users/robert/matlab/fieldtrip-fieldtrip-851478d
  ftver = ftpath(end-6:end);
  
elseif isequal(regexp(ftpath, ['.*\' filesep '[fF]ieldtrip-20[0-9]{6}']), 1)
  % this corresponds with the daily version from the ftp server
  % which results in a ftpath like /Users/robert/matlab/fieldtrip-20160317
  ftver = ftpath(end-7:end);
  
elseif isequal(regexp(ftpath, ['.*\' filesep '[fF]ieldtrip-lite20[0-9]{6}']), 1)
  % this corresponds with the daily version from the ftp server
  % which results in a ftpath like /Users/robert/matlab/fieldtrip-20160317
  ftver = ftpath(end-7:end);
  
else
  % get it from the Contents.m file in the FieldTrip directory
  if ~isdeployed
    tmp = ver(ftpath);
    ftver = tmp.Version;
  else
    ftver = 'deployed';
  end
end % if issvn, isgit or otherwise

if ~nargout
  switch command
    case 'revision'
      fprintf('\nThis is FieldTrip, revision %s.\n\n', ftver);
    case 'branch'
      fprintf('\nThis is FieldTrip, branch %s.\n\n', branch);
    case 'clean'
      fprintf('\nThis is FieldTrip, clean %s.\n\n', clean);
  end
else
  switch command
    case 'revision'
      varargout = {ftver, ftpath};
    case 'branch'
      varargout = {branch};
    case 'clean'
      varargout = {clean};
  end
end
