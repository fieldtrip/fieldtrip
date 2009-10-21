function [tri] = projecttri(pnt, method)

% PROJECTTRI makes a closed triangulation of a list of vertices by
% projecting them onto a unit sphere and subsequently by constructing
% a convex hull triangulation.
%
% Use as
%   [tri] = projecttri(pnt, method)
% The optional method argument can be 'convhull' (default) or 'delaunay'.

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: projecttri.m,v $
% Revision 1.3  2009/09/23 09:05:39  roboos
% added delaunay method
%
% Revision 1.2  2007/05/08 07:36:42  roboos
% updated help
%
% Revision 1.1  2006/12/12 11:27:45  roboos
% created subfunction into a seperate function
%

if nargin<2
  method = 'convhull';
end

switch method
  case 'convhull'
    ori = (min(pnt) + max(pnt))./2;
    pnt(:,1) = pnt(:,1) - ori(1);
    pnt(:,2) = pnt(:,2) - ori(2);
    pnt(:,3) = pnt(:,3) - ori(3);
    nrm = sqrt(sum(pnt.^2, 2));
    pnt(:,1) = pnt(:,1)./nrm;
    pnt(:,2) = pnt(:,2)./nrm;
    pnt(:,3) = pnt(:,3)./nrm;
    tri = convhulln(pnt);
  case 'delaunay'
    % make a 2d triangulation of the projected points using delaunay
    prj = elproj(pnt);
    tri = delaunay(prj(:,1), prj(:,2));
  otherwise
    error('unsupported method');
end



