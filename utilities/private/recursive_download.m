function recursive_download(webLocation, localFolder)

% RECURSIVE_DOWNLOAD downloads a complete directory from a RESTful web service
%
% Use as
%   recursive_download(webLocation, localFolder)
%
% See also WEBREAD, WEBSAVE, UNTAR, UNZIP, GUNZIP

% Copyright (C) 2023, Konstantinos Tsilimparis
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

% Read the HTML content of the URL
htmlContent = webread(webLocation);
pattern = '<a href="([^"]+)">';
matches = regexp(htmlContent, pattern, 'tokens');

% Iterate over the matches
for i = 2:numel(matches) % Ignore i=1, which is the parent directory link: '../'
  item = matches{i}{1};

  if endsWith(item, '/') % It is a folder
    % Create the necessary directories if they do not exist
    subfolder = fullfile(localFolder, item);
    if ~isfolder(subfolder)
      mkdir(subfolder);
    end

    % Recursively download the subfolder
    subWebLocation = strcat(webLocation, '/', item);
    recursive_download(subWebLocation, subfolder);

  else % It is a file
    % Create the necessary directories if they do not exist
    if ~isfolder(localFolder)
      mkdir(localFolder);
    end

    % Download the file
    fileUrl = strcat(webLocation, '/', item);
    localFilePath = fullfile(localFolder, item);
    websave(localFilePath, fileUrl);
  end
end
