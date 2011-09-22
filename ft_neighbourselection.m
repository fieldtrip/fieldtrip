function neighbours = ft_neighbourselection(varargin)

% FT_NEIGHBOURSELECTION is obsolete (replaced by FT_PREPARE_NEIGHBOURS).
%
% See also FT_NEIGHBOURPLOT, FT_PREPARE_NEIGHBOURS


% Copyright (C) 2006-2011, Eric Maris, J?rn M. Horschig, Robert Oostenveld
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


warning('ft_neighbourselection is only a compatibility wrapper, which will soon be removed. Please instead call ft_prepare_neighbours.');
neighbours = ft_prepare_neighbours(varargin{:});

