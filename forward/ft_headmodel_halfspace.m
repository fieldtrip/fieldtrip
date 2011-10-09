function vol = ft_headmodel_halfspace(geom, Pc, varargin)

% FT_HEADMODEL_HALFSPACE creates an EEG volume conduction model that
% is described with an infinite conductive halfspace. You can think
% of this as a plane with on one side a infinite mass of conductive
% material (e.g. water) and on the other side non-conductive material
% (e.g. air).
%
% Use as
%    vol = ft_headmodel_halfspace(geom, Pc, ...)
% where
%   geom.pnt = Nx3 vector specifying N points through which a plane is fitted 
%   Pc       = 1x3 vector specifying the spatial position of a point lying in the conductive halfspace 
%              (this determines the plane normal's direction)
%
% Additional optional arguments should be specified as key-value pairs and can include
%   'sourcemodel'  = string, 'monopole' or 'dipole' (default = 'dipole')
%   'conductivity' = number,  conductivity value of the conductive halfspace (default = 1)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

model = ft_getopt(varargin, 'sourcemodel', 'dipole');
cond  = ft_getopt(varargin, 'conductivity'); 

if isempty(cond)
  warning('Conductivity was not specified, using 1');
  cond = 1;
end

% the description of this volume conduction model consists of the
% description of the plane, and a point in the void halfspace

if isstruct(geom) && isfield(geom,'pnt')
  pnt = geom.pnt;
elseif size(geom,2)==3
  pnt = geom;
else
  error('incorrect specification of the geometry');
end

% fit a plane to the points
[N,P] = fit_plane(pnt);

% checks if Pc is in the conductive part. If not, flip
incond = acos(dot(N,(Pc-P)./norm(Pc-P))) > pi/2;
if ~incond
  N = -N;
end

vol       = [];
vol.cond  = cond;
vol.pnt   = P(:)'; % a point that lies on the plane that separates the conductive tissue from the air
vol.ori   = N(:)'; % a unit vector pointing towards the air
vol.ori   = vol.ori/norm(vol.ori);

if strcmpi(model,'dipole')
  vol.type  = 'halfspace';    
elseif strcmpi(model,'monopole')
  vol.type  = 'halfspace_monopole';    
else
  error('unknow method')
end

function [N,P] = fit_plane(X)
% Fits a plane through a number of points in 3D cartesian coordinates
P = mean(X,1);  % the plane is spanned by this point and by a normal vector
X = bsxfun(@minus,X,P);
[u, s, v] = svd(X, 0);
N = v(:,3); % orientation of the plane, can be in either direction


