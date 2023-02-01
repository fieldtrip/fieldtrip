function headmodel = ft_headmodel_bemcp(mesh, varargin)

% FT_HEADMODEL_BEMCP creates a volume conduction model of the head
% using the boundary element method (BEM) for EEG. This function
% takes as input the triangulated surfaces that describe the boundaries
% and returns as output a volume conduction model which can be used
% to compute leadfields.
%
% The implementation of this function is based on Christophe Phillips'
% MATLAB code, hence the name "bemcp".
%
% Use as
%   headmodel = ft_headmodel_bemcp(mesh, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   conductivity     = vector, conductivity of each compartment
%   checkmesh        = 'yes' or 'no'
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2010, Robert Oostenveld, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

ft_hastoolbox('bemcp', 1);

% get the optional input arguments
conductivity    = ft_getopt(varargin, 'conductivity');
checkmesh       = ft_getopt(varargin, 'checkmesh', 'yes');

% convert to Boolean value
checkmesh = istrue(checkmesh);

if isfield(mesh, 'bnd')
  mesh = mesh.bnd;
end

% replace pnt with pos
mesh = fixpos(mesh);

% determine the number of compartments
numboundaries = length(mesh);

% bemcp expects surface normals to point outwards;
% this checks and corrects if needed
for i=1:numboundaries
  switch surface_orientation(mesh(i).pos, mesh(i).tri)
    case 'outward'
      % this is ok
    case 'inward'
      ft_warning('flipping mesh %d', i);
      mesh(i).tri = fliplr(mesh(i).tri);
    case 'otherwise'
      ft_error('incorrect mesh %d', i)
  end
end

% ensure that the vertices and triangles are double precision, otherwise the bemcp mex files will crash
for i=1:numboundaries
  mesh(i).pos = double(mesh(i).pos);
  mesh(i).tri = double(mesh(i).tri);
end

% start with an empty volume conductor
headmodel = [];

% determine the number of compartments
numboundaries = length(mesh);

if numboundaries~=3
  ft_error('this only works for three surfaces');
end

% determine the desired nesting of the compartments
order = surface_nesting(mesh, 'insidefirst');

% rearrange boundaries and conductivities
if numboundaries>1 && ~isequal(order(:)', 1:numboundaries)
  fprintf('reordering the boundaries to: ');
  fprintf('%d ', order);
  fprintf('\n');
  % update the order of the compartments
  mesh = mesh(order);
end

if isempty(conductivity)
  ft_warning('no conductivity specified, using default values')
  % brain/skull/skin
  conductivity = [0.33 0.0042 0.33];
  headmodel.cond = conductivity;
else
  if numel(conductivity)~=numboundaries
    ft_error('each compartment should have a conductivity value');
  end
  headmodel.cond = conductivity(order);
end

% do some sanity checks on the meshes
if checkmesh
  for i=1:numboundaries
    ntri = size(mesh(i).tri,1);
    npos = size(mesh(i).pos,1);
    assert(2*(npos-2)==ntri, 'the number of triangles does not match the number of vertices')
  end

  % check for each of the vertices that it falls inside the next surface
  for i=1:(numboundaries-1)
    for j=1:size(mesh(i).pos,1)
      inside = surface_inside(mesh(i).pos(j,:), mesh(i+1).pos, mesh(i+1).tri);
      assert(inside==true, 'vertex %d of surface %d is outside surface %d', j, i, i+1);
    end
  end

  ft_info('the meshes are closed and properly nested')
end % if checkmesh

headmodel.bnd          = mesh;
headmodel.skin_surface = numboundaries;
headmodel.source       = 1;

if headmodel.skin_surface~=3
  ft_error('the third surface should be the skin');
end

if headmodel.source~=1
  ft_error('the first surface should the source compartment');
end

% just to be sure
clear mesh

% Build Triangle 4th point
headmodel = triangle4pt(headmodel);

% 2. BEM model estimation, only for the scalp surface

defl =[ 0 0 1/size(headmodel.bnd(headmodel.skin_surface).pos,1)];
% ensure deflation for skin surface, i.e. average reference over skin

% NOTE:
% Calculation proceeds by estimating each submatrix C_ij and combine them.
% There are 2 options:
% - calculating the matrices once, as it takes some time, keep them in
%   memory and use them the 2-3 times they're needed.
% - calculating the matrices every time they're needed, i.e. 2-3 times
% The former option is faster but requires more memory space as up to *8*
% square matrices of size C_ij have to be kept in memory at once.
% The latter option requires less memory, but would take much more time to
% estimate.
% This faster but memory hungry solution is implemented here.

% Deal first with surface 1 and 2 (inner and outer skull
%--------------------------------

% NOTE:
% C11st/C22st/C33st are simply the matrix C11/C22/C33 minus the identity
% matrix, i.e. C11st = C11-eye(N)

weight = (headmodel.cond(1)-headmodel.cond(2))/((headmodel.cond(1)+headmodel.cond(2))*2*pi);
C11st  = bem_Cii_lin(headmodel.bnd(1).tri,headmodel.bnd(1).pos, weight,defl(1),headmodel.bnd(1).pnt4);
weight = (headmodel.cond(1)-headmodel.cond(2))/((headmodel.cond(2)+headmodel.cond(3))*2*pi);
C21    = bem_Cij_lin(headmodel.bnd(2).pos,headmodel.bnd(1).pos,headmodel.bnd(1).tri, weight,defl(1));
tmp1   = C21/C11st;

weight = (headmodel.cond(2)-headmodel.cond(3))/((headmodel.cond(1)+headmodel.cond(2))*2*pi);
C12    = bem_Cij_lin(headmodel.bnd(1).pos,headmodel.bnd(2).pos,headmodel.bnd(2).tri, weight,defl(2));
weight = (headmodel.cond(2)-headmodel.cond(3))/((headmodel.cond(2)+headmodel.cond(3))*2*pi);
C22st  = bem_Cii_lin(headmodel.bnd(2).tri,headmodel.bnd(2).pos, weight,defl(2),headmodel.bnd(2).pnt4);
tmp2   = C12/C22st;

% Try to spare some memory:
tmp10 = - tmp2 * C21 + C11st;
clear C21 C11st
tmp11 = - tmp1 * C12 + C22st;
clear C12 C22st

% Combine with the effect of surface 3 (scalp) on the first 2
%------------------------------------------------------------
weight = (headmodel.cond(1)-headmodel.cond(2))/(headmodel.cond(3)*2*pi);
C31    = bem_Cij_lin(headmodel.bnd(3).pos,headmodel.bnd(1).pos,headmodel.bnd(1).tri, weight,defl(1));
%   tmp4   = C31/(- tmp2 * C21 + C11st );
%   clear C31 C21 C11st
tmp4 = C31/tmp10;
clear C31 tmp10

weight = (headmodel.cond(2)-headmodel.cond(3))/(headmodel.cond(3)*2*pi);
C32    = bem_Cij_lin(headmodel.bnd(3).pos,headmodel.bnd(2).pos,headmodel.bnd(2).tri, weight,defl(2));
%   tmp3   = C32/(- tmp1 * C12 + C22st );
%   clear  C12 C22st C32
tmp3 = C32/tmp11;
clear C32 tmp11

tmp5 = tmp3*tmp1-tmp4;
tmp6 = tmp4*tmp2-tmp3;
clear tmp1 tmp2 tmp3 tmp4

% Finally include effect of surface 3 on the others
%--------------------------------------------------
% As the gama1 intermediate matrix is built as the sum of 3 matrices, I can
% spare some memory by building them one at a time, and summing directly
weight = headmodel.cond(3)/((headmodel.cond(1)+headmodel.cond(2))*2*pi);
Ci3    = bem_Cij_lin(headmodel.bnd(1).pos,headmodel.bnd(3).pos,headmodel.bnd(3).tri, weight,defl(3));
gama1  = - tmp5*Ci3; % gama1 = - tmp5*C13;

weight = headmodel.cond(3)/((headmodel.cond(2)+headmodel.cond(3))*2*pi);
Ci3    = bem_Cij_lin(headmodel.bnd(2).pos,headmodel.bnd(3).pos,headmodel.bnd(3).tri, weight,defl(3));
gama1  = gama1 - tmp6*Ci3; % gama1 = - tmp5*C13 - tmp6*C23;

weight = 1/(2*pi);
Ci3    = bem_Cii_lin(headmodel.bnd(3).tri,headmodel.bnd(3).pos, weight,defl(3),headmodel.bnd(3).pnt4);
gama1  = gama1 - Ci3; % gama1 = - tmp5*C13 - tmp6*C23 - C33st;
clear Ci3

% Build system matrix
%--------------------
i_gama1 = inv(gama1);
headmodel.mat = [i_gama1*tmp5 i_gama1*tmp6 i_gama1];

% remember the type
headmodel.type = 'bemcp';
