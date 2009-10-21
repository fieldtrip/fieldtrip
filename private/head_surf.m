function [pnt,tri]=head_surf(vol,grad,flag);

% HEAD_SURF determines the head surface from a multisphere headmodel
% and a set of coils in the MEG gradiometer system and returns a
% triangulation
% 
% Use as
%   [pnt, tri] = head_surf(vol, grad, flag) 
% where
%   grad	gradiometer definition
%   vol		multisphere volume conductor definition
%
% If flag=1 the lower rim of the helmet-shaped head surface will
% be shifted downward, if flag=0 it will not be shifted downward.
% The default is flag=1.
%
% The head surface triangulation can be used in combination with
% find_inside_vol to determine which voxels of a beamformer scan 
% are located within the head (region of interest).

% Copyright (C) Jan-Matthijs Schoffelen
%
% $Log: head_surf.m,v $
% Revision 1.7  2005/11/01 09:53:17  roboos
% added optional input argument 'flag', which determines whether the lower rim
% will be shifted downward
%
% Revision 1.6  2005/09/29 00:37:33  roboos
% made the code independent of the gradiometer and volume label
% changed the shift of the lower rim, now with distance that is 1/4 of the radius
% updated the help
%
% Revision 1.5  2004/01/26 09:03:12  roberto
% replaced complex 2x2D delaunay with single 3D convex hull triangulation
% better coverage of inferio-frontal region
%
% Revision 1.4  2003/11/03 14:43:25  roberto
% included search for matching labels in volume and gradiometer
%
% Revision 1.3  2003/10/07 10:43:02  roberto
% head surface is now only computed from M channels and not reference channels
%
% Revision 1.2  2003/06/03 08:30:32  roberto
% *** empty log message ***
%
% Revision 1.1  2003/04/17 14:57:32  roberto
%

if nargin<3
  flag = 1;
end

Nchans = size(grad.tra, 1);
Ncoils = size(grad.tra, 2);

% for each coil, determine a surface point using the corresponding sphere
vec = grad.pnt - vol.o;
nrm = sqrt(sum(vec.^2,2));
vec = vec ./ [nrm nrm nrm];
pnt = vol.o + vec .* [vol.r vol.r vol.r];

%  make a triangularization that is needed to find the rim of the helmet
prj = elproj(pnt);
tri = delaunay(prj(:,1),prj(:,2));

% find the lower rim of the helmet shape
[pnt,line] = find_mesh_edge(pnt, tri);
edgeind    = unique(line(:));

if flag
  % shift the lower rim of the helmet shape down with approximately 1/4th of its radius
  dist = mean(sqrt(sum((pnt - repmat(mean(pnt,1), Ncoils, 1)).^2, 2)));
  dist = dist/4;
  pnt(edgeind,3) = pnt(edgeind,3) - dist; 
end

% use matlab triangulation algorithm to determine convex hull, this is the
% final triangulation which makes a nice headshape
tri = convhulln(pnt);
