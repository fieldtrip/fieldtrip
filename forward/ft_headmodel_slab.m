function headmodel = ft_headmodel_slab(mesh1, mesh2, Pc, varargin)

% FT_HEADMODEL_SLAB creates an EEG volume conduction model that
% is described with an infinite conductive slab. You can think
% of this as two parallel planes containing a mass of conductive
% material (e.g. water) and externally to them a non-conductive material
% (e.g. air).
%
% Use as
%   headmodel = ft_headmodel_slab(mesh1, mesh2, Pc, varargin)
% where
%   mesh1.pos = Nx3 vector specifying N points through which the 'upper' plane is fitted 
%   mesh2.pos = Nx3 vector specifying N points through which the 'lower' plane is fitted 
%   Pc        = 1x3 vector specifying the spatial position of a point lying in the conductive slab 
%              (this determines the plane's normal's direction)
% 
% Optional arguments should be specified in key-value pairs and can include
%   'sourcemodel'  = 'monopole' 
%   'conductivity' = number ,  conductivity value of the conductive halfspace (default = 1)
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

model = ft_getopt(varargin, 'sourcemodel', 'monopole');
cond  = ft_getopt(varargin, 'conductivity'); 

if isempty(cond)
  ft_warning('Conductivity was not specified, using 1');
  cond = 1;
end

% the description of this volume conduction model consists of the
% description of the plane, and a point in the void halfspace

% replace pnt with pos
mesh1 = fixpos(mesh1);
mesh2 = fixpos(mesh2);

if isstruct(mesh1) && isfield(mesh1,'pos')
  pos1 = mesh1.pos;
  pos2 = mesh2.pos;
elseif size(mesh1,2)==3
  pos1 = mesh1;
  pos2 = mesh2;
else
  ft_error('incorrect specification of the geometry');
end

% fit a plane to the points
[N1,P1] = fit_plane(pos1);
[N2,P2] = fit_plane(pos2);

% checks if Pc is in the conductive part. If not, flip
incond = acos(dot(N1,(Pc-P1)./norm(Pc-P1))) > pi/2;
if ~incond
  N1 = -N1;
end
incond = acos(dot(N2,(Pc-P2)./norm(Pc-P2))) > pi/2;
if ~incond
  N2 = -N2;
end

headmodel       = [];
headmodel.cond  = cond;
headmodel.pos1   = P1(:)'; % a point that lies on the plane that separates the conductive tissue from the air
headmodel.ori1   = N1(:)'; % a unit vector pointing towards the air
headmodel.ori1   = headmodel.ori1/norm(headmodel.ori1);
headmodel.pos2   = P2(:)'; 
headmodel.ori2   = N2(:)'; 
headmodel.ori2   = headmodel.ori2/norm(headmodel.ori2);

if strcmpi(model,'monopole')
  headmodel.type  = 'slab_monopole';    
else
  ft_error('unknow method')
end

function [N,P] = fit_plane(X)
% Fits a plane through a number of points in 3D cartesian coordinates
P = mean(X,1);  % the plane is spanned by this point and by a normal vector
X = bsxfun(@minus,X,P);
[u, s, v] = svd(X, 0);
N = v(:,3); % orientation of the plane, can be in either direction

