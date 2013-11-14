function img=s2v(node,face,div)
%
% img=s2v(node,face,div)
%
% shortcut for surf2vol, coverting a surface to a volumetric image
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%	 node: node list of the triangular surface, 3 columns for x/y/z
%	 face: triangle node indices, each row is a triangle
%	 div:  division number along the shortest edge of the mesh (resolution)
%              if not given, div=50
%
% output:
%	 img: a volumetric binary image at position of ndgrid(xi,yi,zi)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

p0=min(node);
p1=max(node);

if(size(node,1)==0 | size(face,1)==0)
  error('node and face can not be empty');
end

if(nargin<3)
  div=50;
end

if(div==0)
  error('div can not be 0');
end

dx=min(p1-p0)/div;

if(dx<=eps)
  error('the input mesh is in a plane');
end

img=surf2vol(node,face,p0(1)-dx:dx:p1(1)+dx,p0(2)-dx:dx:p1(2)+dx,p0(3)-dx:dx:p1(3)+dx);
