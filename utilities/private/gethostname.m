function host = gethostname()

% HOSTNAME returns the hostname of this computer
%
% Use as
%   str = hostname;
%
% See also GETUSERNAME, GETADDRESS

% Copyright (C) 2011, Eelke Spaak
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
persistent previous_argout
if ~isempty(previous_argout)
  host = previous_argout;
  return
end

if (ispc())
  host = getenv('ComputerName');
elseif (isunix())
  host = getenv('HOSTNAME');
  % the HOSTNAME variable is not guaranteed to be set
  if isempty(host)
    [status, host] = system('hostname -s');
    % remove the ^n from the end of the line
    host = deblank(host);
  end
end

host = strtok(host, '.'); % dots in filenames are not allowed by matlab

if (isempty(host))
  host = 'unknownhost';
end

% remember for subsequent calls
previous_argout = host;

