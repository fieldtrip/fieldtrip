function [proj, dist] = plinproj(l1, l2, r, flag)

% PLINPROJ projects a point onto a line or linepiece
%
% Use as
%   [proj, dist] = plinproj(l1, l2, r, flag)
% where l1 and l2 are the begin and endpoint of the linepiece, and r is
% the point that is projected onto the line
%
% the optional flag can be:
%   0 (default)  project the point anywhere on the complete line
%   1            project the point within or on the edge of the linepiece

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
