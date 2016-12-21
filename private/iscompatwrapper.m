function b = iscompatwrapper(funcname)
%ISCOMPATWRAPPER Checks whether the specified function name will invoke a
% compatibility wrapper or not.
%
% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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
% $Id

ftPath = which('ft_defaults');
ftPath = ftPath(1:end-numel('ft_defaults.m'));

x = which(funcname);
b = strcmp(x(1:end-numel(funcname)-2), [ftPath 'compat' filesep]);

end

