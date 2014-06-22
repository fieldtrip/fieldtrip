function [optarg] = ft_cfg2keyval(cfg)

% CFG2KEYVAL converts between a structure and a cell-array with key-value
% pairs which can be used for optional input arguments.
%
% Use as
%   [optarg] = ft_cfg2keyval(cfg)

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

if ~isempty(cfg)
  optarg = [fieldnames(cfg) struct2cell(cfg)]';
  optarg = optarg(:)';
else
  optarg = {};
end
