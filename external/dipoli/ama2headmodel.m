function [headmodel] = ama2headmodel(ama)
  
% AMA2HEADMODEL converts a dipoli structure with boundary geometries
% and a boundary element method transfer matrix to a volume conduction
% model.
%
% Use as
%   headmodel = ama2headmodel(ama)

% Copyright (C) 2008, Robert Oostenveld
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

headmodel = [];
ngeo = length(ama.geo);
for i=1:ngeo
  headmodel.bnd(i).pos = ama.geo(i).pos;
  headmodel.bnd(i).tri = ama.geo(i).tri;
  headmodel.cond(i) = ama.geo(i).sigmam;
end
headmodel.mat = ama.bi;
npos = size(headmodel.mat,2);
if size(headmodel.mat,1)<npos
  headmodel.mat(npos, npos) = 0; % it should be a square matrix
end
headmodel.mat  = headmodel.mat;
headmodel.type = 'dipoli';
