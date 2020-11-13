function [event]=ft_fetch_event(data)

% FT_FETCH_EVENT mimics the behavior of FT_READ_EVENT, but for a FieldTrip
% raw data structure instead of a file on disk.
%
% Use as
%   event = ft_fetch_event(data)
%
% See also FT_READ_EVENT, FT_FETCH_HEADER, FT_FETCH_DATA

% Copyright (C) 2008, Esther Meeuwissen
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

% check whether input is data
data = ft_checkdata(data, 'datatype', 'raw');

% locate the event structure
if (isfield(data, 'cfg'))
  event = ft_findcfg(data.cfg, 'event');
  if ~isstruct(event)
    % this happens if the cfg has been cleaned up with checksize in ft_checkconfig, it can then be 'empty - removed by ft_checkconfig'
    event = [];
  end
else
  event = [];
end
