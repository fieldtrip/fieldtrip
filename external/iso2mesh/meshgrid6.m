function [node,elem]=meshgrid6(varargin)
%
% [node,elem]=meshgrid6(v1,v2,v3,...)
%
% mesh an ND rectangular lattice by splitting 
% each hypercube into 6 tetrahedra
%
% author: John D'Errico
% URL: http://www.mathworks.com/matlabcentral/newsreader/view_thread/107191
% modified by Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
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
%     [node,elem]=meshgrid6(0:5,0:6,0:4);
%     plotmesh(node,elem);
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

% dimension of the lattice
n = length(varargin);

% create a single n-d hypercube
% list of node of the cube itself
vhc=('1'==dec2bin(0:(2^n-1)));
% permutations of the integers 1:n
p=perms(1:n);
nt=factorial(n);
thc=zeros(nt,n+1);
for i=1:nt
  thc(i,:)=find(all(diff(vhc(:,p(i,:)),[],2)>=0,2))';
end

% build the complete lattice
nodecount = cellfun('length',varargin);
if any(nodecount<2)
  error 'Each dimension must be of size 2 or more.'
end
node = lattice(varargin{:});

% unrolled index into each hyper-rectangle in the lattice
ind = cell(1,n);
for i=1:n
ind{i} = 0:(nodecount(i)-2);
end
ind = lattice(ind{:});
k = cumprod([1,nodecount(1:(end-1))]);
ind = 1+ind*k';
nind = length(ind);

offset=vhc*k';
elem=zeros(nt*nind,n+1);
L=(1:nind)';
for i=1:nt
  elem(L,:)=repmat(ind,1,n+1)+repmat(offset(thc(i,:))',nind,1);
L=L+nind;
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
