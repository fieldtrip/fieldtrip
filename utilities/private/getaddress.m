function address = getip(hostname)

% GETADDRESS returns the IP address
%
% Use as
%   address = getaddress();
% or
%   address = getaddress(hostname);
%
% See also GETUSERNAME, GETHOSTNAME

% Copyright (C) 2015, Robert Oostenveld
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

% this is to speed up subsequent calls
persistent previous_argin previous_argout

if nargin==0
  hostname = [];
end

if ~isempty(previous_argout) && isequal(hostname, previous_argin)
  address = previous_argout;
  return
end

if ~isempty(hostname)
  address = java.net.InetAddress.getByName(hostname);
else
  address = java.net.InetAddress.getLocalHost;
end
address = char(address.getHostAddress);

% remember for subsequent calls
previous_argin  = hostname;
previous_argout = address;
