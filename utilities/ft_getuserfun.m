function func = ft_getuserfun(func, prefix)
% FT_GETUSERFUN will search the Matlab path for a function with the
% appropriate name, and return its name. Considered are, in this order:
%  - the name itself, i.e. you get exactly the same func back as you put in;
%  - the name with the specified prefix;
%  - the name with 'ft_' and the specified prefix.
%
% For example, calling FT_GETUSERFUN('general', 'trialfun') might return a
% function named 'general', 'trialfun_general', or 'ft_trialfun_general',
% whichever of those is found first and is not a compatibility wrapper.
%
% func can be a function handle, in which case it is returned as-is.
%
% The returned function is guaranteed to be a callable Matlab function, and
% is not a FieldTrip compatibility wrapper. If no such function is found,
% the empty array [] will be returned.

% Copyright (C) 2012, Eelke Spaak
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

if isa(func, 'function_handle') || isfunction(func)
  % treat function (or function name) as-is, this is the default
elseif isfunction([prefix '_' func]) && ~iscompatwrapper([prefix '_' func])
  func = [prefix '_' func];
elseif isfunction(['ft_' prefix '_' func])
  func = ['ft_' prefix '_' func];
else
  warning(['the specified function ''' func ''' could not be found']);
  func = [];
end