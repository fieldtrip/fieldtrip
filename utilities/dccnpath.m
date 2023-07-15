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

% alternative1 allows to test with local files in the present working directory
% this is often convenient when initially setting up a new test script while the data is not yet uploaded
[p, f, x] = fileparts(alternative0);
alternative1 = [f x];

if strcmp(alternative1, 'test')
  % this should not be used when the filename is only "test", since that is too generic
  alternative1 = '';
elseif endsWith(pwd, alternative1)
  % exist(alternative1, 'file') also returns 7 in case the present working directory name matches alternative1
  % this is not the match we want
  alternative1 = '';
end

% alternative2 allows for the user to specify the path in the global ft_default variable
% see https://github.com/fieldtrip/fieldtrip/issues/1998
if isfield(ft_default, 'dccnpath')
    if ~ispc
        alternative2 = strrep(alternative0, '/home/common/matlab/fieldtrip', ft_default.dccnpath);
    else
        alternative2 = strrep(alternative0, 'H:\common\matlab\fieldtrip', ft_default.dccnpath);
    end
else
  alternative2 = '';
end

if exist(alternative1, 'file')
  warning('using present working directory %s', alternative1);
  filename = alternative1;
elseif isfield(ft_default, 'dccnpath')
  warning('using local copy %s ', alternative2);

  % Check if filename exists. If not, download it (Todo: after permission)    
  if ~exist(alternative2, 'file')
 
     % Public data are downloaded from https://download.fieldtriptoolbox.org
     if contains(alternative2, 'ftp')
        % Find the right path to the HTTPS download server
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
    else
    % Todo: FTP connection of private data
    end          
  end
  filename = alternative2;
else
  warning('using default DCCN path %s', alternative0);
  filename = alternative0;
end
end
