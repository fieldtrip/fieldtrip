function vol = ft_headmodel_strip(geom1, geom2, Pc, varargin)

% FT_HEADMODEL_STRIP creates an EEG volume conduction model that
% is described with an infinite conductive strip. You can think
% of this as two parallel planes containing a mass of conductive
% material (e.g. water) and externally to them a non-conductive material
% (e.g. air).
%   geom1.pnt = Nx3 vector specifying N points through which the 'upper' plane is fitted 
%   D         = distance of the parallel planes
%   Pc        = 1x3 vector specifying the spatial position of a point lying in the conductive strip 
%              (this determines the plane's normal's direction)
% 
%   Additional optional arguments include:
%   'sourcemodel'  = 'monopole' or 'dipole' (default)
%   'conductivity' = number ,  conductivity value of the conductive halfspace (default = 1)
% 
% Use as
%   vol = ft_headmodel_strip(geom1, geom2, Pc, varargin)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

model = keyval('sourcemodel',  varargin); if isempty(model), model='monopole'; end
cond  = keyval('conductivity', varargin); 
if isempty(cond), cond = 1; warning('Unknown conductivity value (set to 1)'); end

% the description of this volume conduction model consists of the
% description of the plane, and a point in the void halfspace

if isstruct(geom1) && isfield(geom1,'pnt')
  pnt1 = geom1.pnt;
  pnt2 = geom2.pnt;
elseif size(geom1,2)==3
  pnt1 = geom1;
  pnt2 = geom2;
else
  error('incorrect specification of the geometry');
end

% fit a plane to the points
[N1,P1] = fit_plane(pnt1);
[N2,P2] = fit_plane(pnt2);

% checks if Pc is in the conductive part. If not, flip
incond = acos(dot(N1,(Pc-P1)./norm(Pc-P1))) > pi/2;
if ~incond
  N1 = -N1;
end
incond = acos(dot(N2,(Pc-P2)./norm(Pc-P2))) > pi/2;
if ~incond
  N2 = -N2;
end

vol       = [];
vol.cond  = cond;
vol.pnt1   = P1(:)'; % a point that lies on the plane that separates the conductive tissue from the air
vol.ori1   = N1(:)'; % a unit vector pointing towards the air
vol.ori1   = vol.ori1/norm(vol.ori1);
vol.pnt2   = P2(:)'; 
vol.ori2   = N2(:)'; 
vol.ori2   = vol.ori2/norm(vol.ori2);

if strcmpi(model,'monopole')
  vol.type  = 'strip_monopole';    
else
  error('unknow method')
end

function [N,P] = fit_plane(X)
% Fits a plane through a number of points in 3D cartesian coordinates
P = mean(X,1);  % the plane is spanned by this point and by a normal vector
X = bsxfun(@minus,X,P);
[u, s, v] = svd(X, 0);
N = v(:,3); % orientation of the plane, can be in either direction

