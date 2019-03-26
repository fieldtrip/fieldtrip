function M = prepare_mesh_fittemplate(headshape,template)

% PREPARE_MESH_FITTEMPLATE computes an affine transformation matrix between 2 point clouds 
%
% This function relies on cpd toolbox from  Myronenko, see https://sites.google.com/site/myronenko/research/cpd
%
%
% See also FT_PREPARE_MESH

% Copyright (C) 2019, Simon Homoelle
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

% add toolbox cpd
ft_hastoolbox('cpd', 2);

%
opt.corresp = 0;
opt.method  = 'affine';
opt.max_it = 100;
opt.fgt = 0;
opt.tol = 10e-12;
opt.outliers = 0.0;
opt.outliers = 0;
[transform,~] = cpd_register(headshape,template, opt);

M = eye(4,4);
M(1:3,1:3) = transform.R;
M(1:3,4)   = transform.t;

end
