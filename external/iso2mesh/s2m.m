function [node, elem, face] = s2m(v, f, keepratio, maxvol, method, regions, holes, varargin)
%
% [node,elem,face]=s2m(v,f,keepratio,maxvol,method)
% [node,elem,face]=s2m(v,f,keepratio,maxvol,'tetgen',regions,holes)
%
% volumetric mesh generation from a closed surface, shortcut for surf2mesh
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% inputs and outputs are similar to those defined in surf2mesh
%
% if method='cgalpoly', s2m will call cgals2m and keepratio should be a
% structure (as the 'opt' input in cgals2m)
%
% input default values:
%       method: if ignored, iso2mesh uses surf2mesh ('tetgen') to do the
%               tetrahedral mesh generation
%       regions,holes: if ignored, iso2mesh assumes both are empty
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if (nargin >= 5 && strcmp(method, 'cgalpoly'))
    [node, elem, face] = cgals2m(v, f, keepratio, maxvol);
    return
end
if (nargin <= 5)
    regions = [];
end
if (nargin <= 6)
    holes = [];
end
if (nargin >= 5)
    if (nargin >= 8)
        [node, elem, face] = surf2mesh(v, f, [], [], keepratio, maxvol, regions, holes, 0, method, varargin{:});
    else
        [node, elem, face] = surf2mesh(v, f, [], [], keepratio, maxvol, regions, holes, 0, method);
    end
else
    [node, elem, face] = surf2mesh(v, f, [], [], keepratio, maxvol, regions, holes);
end
