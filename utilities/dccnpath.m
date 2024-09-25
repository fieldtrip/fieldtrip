function filename = dccnpath(filename)

% DCCNPATH manages the filename and path for test files. It helps to locate and read
% test file from Linux, Windows or macOS computers both inside and outside the DCCN.
%
% Use as
%   filename = dccnpath(filename)
% where the input filename corresponds to the test data on the DCCN cluster and the
% output filename corresponds to the local file including the full path where the
% test data is available.
%
% The test data location on the DCCN cluster is '/home/common/matlab/fieldtrip/data'
% and the specification of the input filename MUST start with this.
%
% This function will search-and-replace the location on the DCCN cluster by the
% location that applies to your computer. If needed, it will replace '/home' by 'H:'
% and will replace forward by backward slashes.
%
% In case you have a local copy of the data, or if you are inside the DCCN and have
% mounted the '/home' drive on another letter than 'H:', you should override the
% default location using
%    global ft_default
%    ft_default.dccnpath = '/your/copy';
%
% If you DO HAVE a local copy, it should contain a directory with the name 'ftp'. The 
% content of the ftp directory should match that on the FieldTrip download server, 
% for example '/your/copy/ftp/test/ctf'.
%
% If you DO NOT have a local copy and do not define ft_default.dccnpath manually,
% then this function will automatically try to download the publicly available data 
% to a temporary directory.
%
% See also WHICH, WEBSAVE
% Copyright (C) 2012-2024, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

global ft_default

% it must always start with this
assert(startsWith(filename, '/home/common/matlab/fieldtrip/data'));

% alternative0 applies inside the DCCN network when using the standard data location
if ispc
  alternative0 = strrep(filename, '/home', 'H:');
  alternative0 = strrep(alternative0, '/', '\');
else
  alternative0 = strrep(filename, 'H:', '/home');
  alternative0 = strrep(alternative0, '\', '/');
end

if exist(alternative0, 'file') || exist(alternative0, 'dir')
  ft_notice('using default DCCN path %s', alternative0);
  filename = alternative0;
  return
end

% the simple alternative1 applies with a local file in the present working directory
% this is often convenient when initially setting up a new test script while the data is not yet uploaded
[p, f, x] = fileparts(alternative0);
alternative1 = [f x];

skip = {
  % subdirectories like fieldtrip/xxx
  'bin'
  'compat'
  'connectivity'
  'contrib'
  'external'
  'fileio'
  'forward'
  'inverse'
  'plotting'
  'preproc'
  'private'
  'qsub'
  'realtime'
  'specest'
  'src'
  'statfun'
  'template'
  'test'
  'trialfun'
  'utilities'
  % subdirectories like fieldtrip/template/xxx
  'dewar'
  'headmodel'
  'sourcemodel'
  'anatomy'
  'electrode'
  'layout'
  'atlas'
  'gradiometer'
  'neighbours'
  };
if any(strcmp(alternative1, skip))
  % this should not be used for subdirectories underneath fieldtrip
  warning('the simple alternative1 cannot be used for "%s"', alternative1)
  alternative1 = '';
end

if exist(alternative1, 'file') || exist(alternative1, 'dir')
  if ~isempty(x) && ~isequal(x, '.ds')   % if alternative1 is a file
    filenamepath = which(alternative1); % also output the path that alternative1 is located at
    ft_notice('using present working directory %s', filenamepath);
    filename = filenamepath;
    return;
  else                 % if alternative1 is a folder
    searchPath = path; % get the MATLAB search path
    directories = strsplit(searchPath, pathsep); % split the search path into individual directories

    folderPath = '';
    for i = 1:numel(directories)
      currentFolder = directories{i};
      if contains(currentFolder, alternative1)
        folderPath = currentFolder;
        break
      end
    end
    warning('using present working directory %s', folderPath);
    filename = folderPath;
    return
  end
end

% alternative2 applies when ft_default.dccnpath is specified
% see https://github.com/fieldtrip/fieldtrip/issues/1998

if ~isfield(ft_default, 'dccnpath') || isempty(ft_default.dccnpath)
  ft_notice('using a temporary directory')
  ft_default.dccnpath = tempdir;
end

% we do not want it to end with a '/' or '\'
ft_default.dccnpath = strip(ft_default.dccnpath, 'right', '/');
ft_default.dccnpath = strip(ft_default.dccnpath, 'right', '\');

% alternative0 is the same as the input filename, but potentially updated for windows
if ~ispc
  alternative2 = strrep(alternative0, '/home/common/matlab/fieldtrip/data', ft_default.dccnpath);
else
  alternative2 = strrep(alternative0, 'H:\common\matlab\fieldtrip\data', ft_default.dccnpath);
end

if exist(alternative2, 'file')
  ft_notice('using local copy %s ', alternative2);
  filename = alternative2;
  return

elseif isfolder(alternative2) && ~isemptydir(alternative2)
  ft_notice('using local copy %s ', alternative2);
  filename = alternative2;
  return

else
  % if the file doesn't exist or the folder is empty, then download test data
  % see also UNTAR, UNZIP, GUNZIP, which can download on the fly

  if contains(alternative0, 'data/test') || contains(alternative0, 'data\test')
    error('the test data are private and can not be downloaded from https://download.fieldtriptoolbox.org')
  end

  % public data are downloaded from https://download.fieldtriptoolbox.org
  % so, we need to find the right path to the HTTPS download server
  pattern = 'ftp(.*)';
  datadir = regexp(alternative0, pattern, 'tokens', 'once');
  if iscell(datadir)
    datadir = datadir{1};
  end
  weblocation = strcat('https://download.fieldtriptoolbox.org', datadir);
  weblocation = strrep(weblocation, '\', '/');

  urlContent = webread(weblocation, weboptions('ContentType', 'text'));

  if contains(urlContent, '<html')
    % the URL corresponds to a folder
    ft_notice('downloading recursively from %s', weblocation);
    recursive_download(weblocation, alternative2);
    filename = alternative2;

  else
    % the URL corresponds to a file
    % create the necessary directories if they do not exist
    [p, f, x] = fileparts(alternative2);
    if ~isfolder(p)
      mkdir(p);
    end

    if isfolder(alternative2)
      [p, f, x] = fileparts(alternative0);
      alternative2 = fullfile(alternative2, [f x]);
    end

    ft_notice('downloading recursively from %s', weblocation);
    websave(alternative2, weblocation);
    filename = alternative2;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = isemptydir(dirName)
dirContents = dir(dirName);
bool = ~isempty(dirContents(~ismember({dirContents.name}, {'.', '..'})));
