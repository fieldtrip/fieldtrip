function [node,face,elem]=meshcylinders(c0, v, len, varargin)
%
% [node,face]=meshcylinders(c0, v, len, r,tsize,maxvol,ndiv)
%    or
% [node,face,elem]=meshacylinder(c0, v, len, r, tsize,maxvol,ndiv)
% [nplc,fplc]=meshacylinder(c0, v, len,r,0,0,ndiv);
%
% create the surface and (optionally) tetrahedral mesh of a 3D cylinder
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input: 
%   c0, cylinder list axis's starting point
%   v: directional vector of the cylinder
%   len: a scalar or a vector denoting the length of each 
%        cylinder segment along the direction of v
%   tsize, maxvol, ndiv: please see the help for meshacylinder for details
%
% output:
%   node, face, elem: please see the help for meshacylinder for details
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

len=cumsum(len);
[ncyl,fcyl]=meshacylinder(c0,c0+v*len(1),varargin{:});

for i=2:length(len)
   [ncyl1,fcyl1]=meshacylinder(c0+v*len(i-1),c0+v*len(i),varargin{:});
   fcyl1=cellfun(@(x) {x{1}+size(ncyl,1),x{2}}, fcyl1, 'UniformOutput', false);
   ncyl=[ncyl; ncyl1];
   if(i==1)
       fcyl1=fcyl1(1:end-1);
   else
       fcyl1={fcyl1{1:end-2},fcyl1{end}};
   end
   fcyl={fcyl{:}, fcyl1{:}};
end

[ncyl,I,J]=unique(round(ncyl*1e10),'rows');
ncyl=ncyl*1e-10;
fcyl=cellfun(@(x) {J(x{1})',x{2}}, fcyl, 'UniformOutput', false);

tsize=varargin{2};
maxvol=varargin{3};

if(nargout==2 && tsize==0.0 && maxvol==0.0)
    node=ncyl;
    face=fcyl;
    return;
end
if(nargin==3)
    tsize=len/10;
end
if(nargin<5)
    maxvol=tsize*tsize*tsize;
end

centroid=cumsum([0 len(1:end-1)])+len/2;   % define the centroids of each cylinder segment
seeds=repmat(c0(:)',length(len),1)+repmat(v(:)',length(len),1).*repmat(centroid(:),1,3);
[node,elem,face]=surf2mesh(no,fc,min(no),max(no),1,maxvol,seeds,[],0);
