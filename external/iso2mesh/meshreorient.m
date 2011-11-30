function newelem=meshreorient(node,elem)
%
% newelem=meshreorient(node,elem)
%
% reorder nodes in a surface or tetrahedral mesh to ensure all
% elements are oriented consistently
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2010/05/05
%
% input:
%    node: list of nodes
%    elem: list of elements (each row are indices of nodes of each element)
%
% output:
%    newelem: the element list with consistent ordering
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

% calculate the canonical volume of the element (can be a 2D or 3D)
vol=elemvolume(node,elem,'signed');

% make sure all elements are positive in volume
idx=find(vol<0);
elem(idx,[end-1,end])=elem(idx,[end,end-1]);
newelem=elem;
