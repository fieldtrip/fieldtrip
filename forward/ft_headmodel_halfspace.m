function headmodel = ft_headmodel_halfspace(mesh, Pc, varargin)

% FT_HEADMODEL_HALFSPACE creates an EEG volume conduction model that
% is described with an infinite conductive halfspace. You can think
% of this as a plane with on one side a infinite mass of conductive
% material (e.g. water) and on the other side non-conductive material
% (e.g. air).
%
% Use as
%    headmodel = ft_headmodel_halfspace(mesh, Pc, ...)
% where
%   mesh.pos = Nx3 vector specifying N points through which a plane is fitted 
%   Pc       = 1x3 vector specifying the spatial position of a point lying in the conductive halfspace 
%              (this determines the plane normal's direction)
%
% Additional optional arguments should be specified as key-value pairs and can include
%   'sourcemodel'  = string, 'monopole' or 'dipole' (default = 'dipole')
%   'conductivity' = number,  conductivity value of the conductive halfspace (default = 1)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

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
% $Id$

model = ft_getopt(varargin, 'sourcemodel', 'dipole');
cond  = ft_getopt(varargin, 'conductivity'); 

if isempty(cond)
  ft_warning('Conductivity was not specified, using 1');
  cond = 1;
end

% the description of this volume conduction model consists of the
% description of the plane, and a point in the void halfspace

if isstruct(mesh) && isfield(mesh,'pos')
  pos = mesh.pos;
elseif size(mesh,2)==3
  pos = mesh;
else
  ft_error('incorrect specification of the geometry');
end

% fit a plane to the points
[N,P] = fit_plane(pos);

% checks if Pc is in the conductive part. If not, flip
incond = acos(dot(N,(Pc-P)./norm(Pc-P))) > pi/2;
if ~incond
  N = -N;
end

headmodel       = [];
headmodel.cond  = cond;
headmodel.pos   = P(:)'; % a point that lies on the plane that separates the conductive tissue from the air
headmodel.ori   = N(:)'; % a unit vector pointing towards the air
headmodel.ori   = headmodel.ori/norm(headmodel.ori);

if strcmpi(model,'dipole')
  headmodel.type  = 'halfspace';    
elseif strcmpi(model,'monopole')
  headmodel.type  = 'halfspace_monopole';    
else
  ft_error('unknow method')
end

function [N,P] = fit_plane(X)
% Fits a plane through a number of points in 3D cartesian coordinates
P = mean(X,1);  % the plane is spanned by this point and by a normal vector
X = bsxfun(@minus,X,P);
[u, s, v] = svd(X, 0);
N = v(:,3); % orientation of the plane, can be in either direction


