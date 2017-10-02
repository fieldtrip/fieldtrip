function [IA, IB] = match_val(A, B)

% MATCH_VAL looks for matching values in two arrays of values
% and returns the indices into both the 1st and 2nd list of the matches.
% They will be ordered according to the first input argument.
%
% Use as
%   [sel1, sel2] = match_str(vallist1, vallist2)
%
% See also MATCH_STR

% Copyright (C) 2017, Robert Oostenveld
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

A = A(:);
B = B(:);

[C,IA,IB] = intersect(A, B, 'rows', 'stable');
