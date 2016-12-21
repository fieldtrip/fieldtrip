function ft_create_buffer(port)

% FT_CREATE_BUFFER starts the thread with the TCP server attached to the local
% MATLAB instance. The TCP server will listen to the specified network
% port, and accept incoming read and write requests.
%
% Use as
%   ft_create_buffer(port)
% where port is the TCP port to which the server listens. The default port 
% number is 1972.
% 
% See also FT_DESTROY_BUFFER

% Copyright (C) 2010, Stefan Klanke
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

if nargin<1
  port = 1972;
end

try
  buffer('tcpserver', 'init', 'localhost', port);
  pause(1);
catch
  if ~isempty(strfind(lasterr, 'thread is already running'))
    warning('thread is already running');
  else
    rethrow(lasterror);
  end
end
