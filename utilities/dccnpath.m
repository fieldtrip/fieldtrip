function filename = dccnpath(filename)

% DCCNPATH manages the filename and path for test files. It helps to locate and read test file from Linux, Windows or macOS computers both inside and outside the DCCN.
% 
% Use as 
%   filename = dccnpath(filename)
% 
% Inputs:
%   filename - The relative path to the test data file within the DCCN cluster.
%              Examples:
%              - For public test data: '/home/common/matlab/fieldtrip/data/ftp/testdata' 
%              - For private test data: '/home/common/matlab/fieldtrip/data/test/testdata' 
% 
% Outputs:
%   filename - string that corresponds to the local path that the test data
%   are saved.
% 
% The default location for FieldTrip and its test data is '/home/common/matlab/fieldtrip'. This function will search-and-replace this string by the location that applies to your computer. It will replace '/home' by 'H:' and will replace forward by backward slashes. 
% 
% In case you have a local copy:
% 
% - You should override the default location by 
%    global ft_default 
%    ft_default.dccnpath = '/your/copy'; 
% 
%    If you do not define ft_default.dccnpath manually then dccnpath will
%    automatically use: 
%    ft_default.dccnpath = tempdir;
% 
% - If the test data are not downloaded inside '/your/copy', then dccnpath
%   will automatically download them.
% 
% Copyright (C) 2012-2023, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL

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


% alternative0 allows to test when in the DCCN cluster
if ispc
  alternative0 = strrep(filename,'/home','H:');
  alternative0 = strrep(alternative0,'/','\');
else
  alternative0 = strrep(filename,'H:','/home');
  alternative0 = strrep(alternative0,'\','/');
end

if exist(alternative0, 'file') || exist(alternative0, 'dir')
    warning('using default DCCN path %s', alternative0);
    filename = alternative0;
    return;
end



% alternative1 allows to test with local files in the present working directory
% this is often convenient when initially setting up a new test script while the data is not yet uploaded
[p, f, x] = fileparts(alternative0);
alternative1 = [f x];

if strcmp(alternative1, 'test')
  % this should not be used when the filename is only "test", since that is too generic
  alternative1 = '';
  warning('You are currently in the test directory, change directory so dccnpath can function properly')
% elseif endsWith(pwd, alternative1)
%   % exist(alternative1, 'file') also returns 7 in case the present working directory name matches alternative1
%   % this is not the match we want
%   alternative1 = '';
elseif strcmp(alternative1, 'ctf')
  % this should not be used when the filename is only "ctf", since that is too generic 
  alternative1 = '';
end

if exist(alternative1, 'file') || exist(alternative1, 'dir') 
    if ~isempty(x) && ~isequal(x,'.ds') % if alternative1 is a file
      filenamepath = which(alternative1); % output the path that alternative1 is located at
      warning('using present working directory %s', filenamepath);
      filename = filenamepath;
      return;
    else  % if alternative1 is a folder
        searchPath = path; % Get the MATLAB search path
        directories = strsplit(searchPath, pathsep); % Split the search path into individual directories
        
        folderPath = '';
        for i = 1:numel(directories)
            currentFolder = directories{i};
            if contains(currentFolder, alternative1)
                folderPath = currentFolder;
                break;
            end
        end
        warning('using present working directory %s', folderPath);
        filename = folderPath;
        return;
    end 
end



% alternative2 allows for the user to specify the path in the global ft_default variable
% see https://github.com/fieldtrip/fieldtrip/issues/1998
if isfield(ft_default, 'dccnpath')    
        alternative2 = ft_default.dccnpath;      
else % save the test data automatically in the tempdir
    ft_default.dccnpath = tempdir;
    ft_default.dccnpath = ft_default.dccnpath(1:end-1); % we do not want the temppdir to end with a '/' or '\'   
    if ~ispc
        alternative2 = strrep(alternative0, '/home/common/matlab/fieldtrip', ft_default.dccnpath);
    else
        alternative2 = strrep(alternative0, 'H:\common\matlab\fieldtrip', ft_default.dccnpath);
    end
end

% if folder or file doesn't exist or is empty, then download test data
dirContents = dir(alternative2);
if ~exist(alternative2, 'file') || ~exist(alternative2, 'dir') || isempty(dirContents(~ismember({dirContents.name}, {'.', '..'}))) 

    if contains(alternative0, 'data/test') || contains(alternative0, 'data\test')
        error('The test data are private and can not be downloaded from the public download server')
    end

   
    % Public data are downloaded from https://download.fieldtriptoolbox.org. So, we need to find the right path to the HTTPS download server
    pattern = 'ftp(.*)';                      
    datadir = regexp(alternative0, pattern, 'tokens', 'once');          
    weblocation = strcat("https://download.fieldtriptoolbox.org", datadir);
    weblocation = strrep(weblocation, '\', '/');
    
    urlContent = webread(weblocation, weboptions('ContentType', 'text')); 
    
    % Check if the URL corresponds to a folder or file
    if contains(urlContent, '<html') % Folder (HTML content) 
          recursiveDownload(weblocation, alternative2);
    else % File
    
          % Create the necessary directories if they do not exist
          [folder, ~, ~] = fileparts(alternative2); 
          if ~isfolder(folder)
              mkdir(folder);
          end
          
          if isfolder(alternative2)
              [~, file, fileextension] = fileparts(alternative0);
              alternative2 = fullfile(alternative2, [file fileextension]);
          end
                         
          websave(alternative2, weblocation);
    end
end 

if exist(alternative2, 'file') || exist(alternative2, 'dir') 
    warning('using local copy %s ', alternative2);
    filename = alternative2;
%     ft_default=rmfield(ft_default,'dccnpath');
    return;
end
    
end

