function [node,elem]=meshgrid5(varargin)
%
% [node,elem]=meshgrid5(v1,v2,v3,...)
%
% mesh an ND rectangular lattice by splitting 
% each hypercube into 5 tetrahedra
%
% author: Qianqian Fang, <q.fang at neu.edu>
% inspired by John D'Errico
% URL: http://www.mathworks.com/matlabcentral/newsreader/view_thread/107191
%
% input:
%    v1,v2,v3,... - numeric vectors defining the lattice in
%                   each dimension.
%                   Each vector must be of length >= 1
%
% output:
%    node - factorial lattice created from (v1,v2,v3,...)
%           Each row of this array is one node in the lattice
%    elem - integer array defining simplexes as references to
%           rows of "node".
%
% example:
%     [node,elem]=meshgrid5(0:5,0:6,0:4);
%     plotmesh(node,elem);
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

% dimension of the lattice
n = length(varargin);
if(n~=3) error('only works for 3D case!'); end

for i=1:n
	v=varargin{i};
	if(mod(length(v),2)==0)
		varargin{i}=linspace(v(1),v(end),length(v)+1);
	end
end

% create a single n-d hypercube
% list of node of the cube itself

cube8=...
[1 4 5 13;1 2 5 11;1 10 11 13;11 13 14 5;11 13 1 5;...
 2 3 5 11;3 5 6 15;15 11 12 3;15 11 14 5;11 15 3 5;...
 4 5 7 13;5 7 8 17;16 17 13 7;13 17 14 5;5 7 17 13;...
 5 6 9 15;5 8 9 17;17 18 15 9;17 15 14 5;17 15 5 9;...
 10 13 11 19;13 11 14 23;22 19 23 13;19 23 20 11;13 11 19 23;...
 11 12 15 21;11 15 14 23;23 21 20 11;23 24 21 15;23 21 11 15;...
 16 13 17 25;13 17 14 23;25 26 23 17;25 22 23 13;13 17 25 23;...
 17 18 15 27;17 15 14 23;26 27 23 17;27 23 24 15;23 27 17 15]';

% build the complete lattice
nodecount = cellfun('length',varargin);
if any(nodecount<2)
	error 'Each dimension must be of size 2 or more.'
end
node = lattice(varargin{:});

[ix,iy,iz]=meshgrid(1:2:nodecount(1)-2,1:2:nodecount(2)-2,1:2:nodecount(3)-2);
ind=sub2ind(nodecount,ix(:),iy(:),iz(:));

nodeshift=[0 1 2 nodecount(1) nodecount(1)+1 nodecount(1)+2 ...
2*nodecount(1) 2*nodecount(1)+1 2*nodecount(1)+2];
nodeshift=[nodeshift,nodeshift+nodecount(1)*nodecount(2),nodeshift+2*nodecount(1)*nodecount(2)];

nc=length(ind);
elem=zeros(nc*40,4);
for i=1:nc
	elem((1:40)+(i-1)*40,:)=reshape(nodeshift(cube8(:)),4,40)'+ind(i);
end

% ======== subfunction ========
function g = lattice(varargin)
% generate a factorial lattice in n variables
n=nargin;
sizes = cellfun('length',varargin);
c=cell(1,n);
[c{1:n}]=ndgrid(varargin{:});
g=zeros(prod(sizes),n);
for i=1:n
	g(:,i)=c{i}(:);
end
