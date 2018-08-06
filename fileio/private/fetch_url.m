function filename = fetch_url(filename)

% FETCH_URL checks the filename and downloads the file to a local copy in
% case it is specified as an Universal Resource Locator. It returns the
% name of the temporary file on the local filesystem.
%
% Use as
%   filename = fetch_url(filename)
%
% In case the filename does not specify an URL, it just returns the original
% filename.

% Copyright (C) 2012 Robert Oostenveld
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

if filetype_check_uri(filename, 'sftp')
  [user, host, filename] = filetype_check_uri(filename);
  [p, f, x] = fileparts(filename);
  p = tempdir;
  try
    mkdir(p);
    cmd = sprintf('sftp %s@%s:%s %s', user, host, filename, fullfile(p, [f x]));
    system(cmd);
    filename = fullfile(p, [f x]);
  end
  % elseif filetype_check_uri(filename, 'http')
  % FIXME the http scheme should be supported using default MATLAB
  % elseif filetype_check_uri(filename, 'ftp')
  % FIXME the http scheme should be supported using default MATLAB
  % elseif filetype_check_uri(filename, 'smb')
  % FIXME the smb scheme can be supported using smbclient
end
