function [proj, dist] = ptriproj(v1, v2, v3, r, flag)

% PTRIPROJ projects a point onto the plane going through a triangle
%
% Use as
%   [proj, dist] = ptriproj(v1, v2, v3, r, flag)
% where v1, v2 and v3 are three vertices of the triangle, and r is 
% the point that is projected onto the plane spanned by the vertices
%
% the optional flag can be:
%   0 (default)  project the point anywhere on the complete plane
%   1            project the point within or on the edge of the triangle

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

