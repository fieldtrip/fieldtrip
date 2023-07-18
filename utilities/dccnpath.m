function filename = dccnpath(filename)

% DCCNPATH manages the filename and path for test files. It helps to locate and read
% test file from Linux, Windows or macOS computers both inside and outside the DCCN.
%
% Use as
%  filename = dccnpath(filename)
%
% The default location for FieldTrip and its test data is '/home/common/matlab/fieldtrip'.
% This function will search-and-replace this string by the location that applies to
% your computer. It will replace '/home' by 'H:' and will replace forward by backward slashes.
%
% In case you have a local copy, you can override the default location by
%   global ft_default
%   ft_default.dccnpath = '/your/copy';
%
% Note that most test scripts expect data located at /home/common/matlab/fieldtrip/data/ftp
% or /home/common/matlab/fieldtrip/data/test, hence you should organize your local
% copy of the data under /your/copy/data/ftp and /your/copy/data/test.
%
% Copyright (C) 2012-2022, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL

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
elseif endsWith(pwd, alternative1)
  % exist(alternative1, 'file') also returns 7 in case the present working directory name matches alternative1
  % this is not the match we want
  alternative1 = '';
end

if exist(alternative1, 'file') || exist(alternative1, 'dir') 
    if ~isempty(x) % if alternative1 is a file
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
            if exist(fullfile(currentFolder, alternative1), 'dir') == 7
                folderPath = fullfile(currentFolder, alternative1);
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
    if ~ispc
        alternative2 = strrep(alternative0, '/home/common/matlab/fieldtrip', ft_default.dccnpath);
    else
        alternative2 = strrep(alternative0, 'H:\common\matlab\fieldtrip', ft_default.dccnpath);
    end
else % save the test data automatically in the tempdir
    ft_default.dccnpath = tempdir;
    ft_default.dccnpath = ft_default.dccnpath(1:end-1); % we do not want the temppdir to end with a '/' or '\'   
    if ~ispc
        alternative2 = strrep(alternative0, '/home/common/matlab/fieldtrip', ft_default.dccnpath);
    else
        alternative2 = strrep(alternative0, 'H:\common\matlab\fieldtrip', ft_default.dccnpath);
    end
end


% if dccnpath's input starts with
% /home/common/matlab/fieldtrip/data/ftp, download test data if it is not already downloaded.
if ~exist(alternative2, 'file') || ~exist(alternative2, 'dir')

    if contains(alternative2, 'data/test') || contains(alternative2, 'data\test')
        warning('Data are private and can not be used if there is no connection to the Donders intranet')
    end

   
    % Public data are downloaded from https://download.fieldtriptoolbox.org. So, we need to find the right path to the HTTPS download server
    pattern = 'ftp(.*)';                      
    datadir = regexp(alternative2, pattern, 'tokens', 'once');          
    weblocation = strcat("https://download.fieldtriptoolbox.org", datadir);
    weblocation = strrep(weblocation, '\', '/');
    
    urlContent = webread(weblocation, weboptions('ContentType', 'text')); % Runs a bit slow for large files
    
    % Check if the URL corresponds to a folder or file
    if contains(urlContent, '<html') % Folder (HTML content) 
          recursiveDownload(weblocation, alternative2);
    else % File
    
          % Create the necessary directories if they do not exist
          [folder, ~, ~] = fileparts(alternative2); 
          if ~isfolder(folder)
              mkdir(folder);
          end
                         
          websave(alternative2, weblocation);
    end
end 

if exist(alternative2, 'file') || exist(alternative2, 'dir') %this "if" just checks if everything went alright. Maybe is not needed
    warning('using local copy %s ', alternative2);
    filename = alternative2;
    return;
end



% if dccnpath=incorrect MATLAB gives automatic error:
% Error using matlab.internal.webservices.HTTPConnector/copyContentToByteArray
% The server returned the status 404 with message "Not Found" in response to the request to URL
% https://download.fieldtriptoolbox.org/tutorial/epilepsy/raw/case1/neuromag/case1_cHPI_ra.fif.

end