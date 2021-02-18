function [output] = volumepad(input, n)

% VOLUMEPAR is a helper function for segmentations. It adds a layer on all sides to
% ensure that the tissue can be meshed all the way up to the edges this also ensures
% that the mesh at the bottom of the neck will be closed.
%
% See also VOLUMEFILLHOLES, VOLUMESMOOTH, VOLUMETHRESHOLD

% Copyrights (C) 2021, Robert Oostenveld
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

if nargin<2 || isempty(n)
  n = 1;
end

dim = size(input);
output = false(dim+2*n);
selx = (1+n):(dim(1)+n);
sely = (1+n):(dim(2)+n);
selz = (1+n):(dim(3)+n);
% insert the original data in the padded boolean volume, the edges remain "false"
output(selx, sely, selz) = input;
