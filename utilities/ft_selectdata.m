function [varargout] = ft_selectdata(varargin)

% FT_SELECTDATA makes a selection in the input data along specific data
% dimensions, such as channels, time, frequency, trials, etc. It can also
% be used to average the data along each of the specific dimensions.
%
% There exist two implementations with a different calling style, i.e. a
% different specification of the input arguments. Eventually all the
% functionality will be merged again. This present FT_SELECTDATA function
% is a wrapper function that will automatically select between the two.
% Please be aware that the old implementation will be phased out.
%
% Use as
%  [data] = ft_selectdata_old(data1, key1, value1, key2, value2, ...)
% or as
%  [data] = ft_selectdata_new(cfg, data, ...)
%
% See also FT_SELECTDATA_OLD, FT_SELECTDATA_NEW

% Copyright (C) 2009-2011, Jan-Mathijs Schoffelen, Robert Oostenveld
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

if nargin==1 || (nargin>2 && ischar(varargin{end-1})) || (isstruct(varargin{1}) && ~ft_datatype(varargin{1}, 'unknown'))
  % this is the OLD calling style, like this
  %   data = ft_selectdata(data, 'key1', value1, 'key2', value2, ...)
  % or with multiple input data structures like this
  %   data = ft_selectdata(data1, data2, 'key1', value1, 'key2', value2, ...)
  [varargout{1:nargout}] = ft_selectdata_old(varargin{:});
else
  % this is the NEW calling style, like this
  %  [data, cfg] = ft_selectdata(cfg, data)
  % or with multiple input data structures like this
  %  [data1, data2, ..., cfg] = ft_selectdata(cfg, data1, data2, ...)
  [varargout{1:nargout}] = ft_selectdata_new(varargin{:});
end
