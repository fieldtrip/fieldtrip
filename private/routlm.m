function [proj] = routlm(v1, v2, v3, la, mu)

% ROUTLM computes the projection of a point from its la/mu parameters
% these equal the "Barycentric" coordinates
%
% Use as
%   [proj] = routlm(v1, v2, v3, la, mu)
% where v1, v2 and v3 are three vertices of the triangle

% Copyright (C) 2002-2009, Robert Oostenveld
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

error('Could not locate the MEX file "%s.%s"', mfilename, mexext);
