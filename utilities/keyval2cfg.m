function [cfg] = keyval2cfg(varargin)

% KEYVAL2CFG converts between a structure and a cell-array with key-value
% pairs which can be used for optional input arguments. 
% 
% Use as
%   [cfg] = keyval2cfg(varargin)

% Copyright (C) 2006, Robert Oostenveld
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

if iscell(varargin) && length(varargin)==1
  varargin = varargin{1};
end

% assign the optional key-value arguments to a configuration structure
var = varargin(1:2:length(varargin));   % get the odd arguments
val = varargin(2:2:length(varargin));   % get the even arguments
cfg = cell2struct(val(:), var(:), 1);
