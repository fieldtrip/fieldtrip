function [node,face,centroids]=latticegrid(varargin)
%
% [node,face,centroids]=latticegrid(xrange,yrange,zrange,...)
%
% generate a 3D lattice
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input: 
%   xrange, yrange, zrange ...: 1D vectors specifying the range of each
%         dimension of the lattice
%
% output:
%   node: the vertices of the 3D lattice
%   face: the list of cell faces of the lattice, including both internal
%         and external facets. By default, face is in the form of a cell
%         array, with each row representing a face. One can use
%         cell2mat(face) to convert it to an array
%   centroids: the centroids of each lattice cell
%
% example:
%    % generate a 3D lattice
%    [node,face,c0]=latticegrid([1 2 4],1:3,1:4);
%    plotmesh(node,face)
%    
%    % mesh the 3D lattice based on the face info
%    [no,el]=surf2mesh(node,face,[],[],1,0.01,c0);
%    figure; plotmesh(no,el)
%
%    % mesh a 2-layer structure using a simple lattice
%    [node,face,c0]=latticegrid([0 10],[0 5],[0 3.5 4]);
%    c0(:,4)=[0.01;0.001];
%    [no,el]=surf2mesh(node,face,[],[],1,[],c0);
%    figure; plotmesh(no,el)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
% 

n=length(varargin);
p=cell(n,1);
[p{:}]=ndgrid(varargin{:});
node=zeros(length(p{1}(:)),n);
for i=1:n
   node(:,i)=p{i}(:);
end
if(nargout==1)
    return;
end

dim=size(p{1});

dd=[dim(1) dim(1)*dim(2)];
onecube=[0 dd(1) dd(1)+1 1; ...
         0 1 dd(2)+1 dd(2); ...
         0 dd(2) dd(2)+dd(1) dd(1)];
onecube=[onecube;onecube+repmat([dd(2);dd(1);1],1,4)];

len=prod(dim(1:3)-1);
face=repmat(onecube,len,1);
[xx,yy,zz]=ndgrid(1:dim(1)-1,1:dim(2)-1,1:dim(3)-1);
idx=sub2ind(dim,xx(:),yy(:),zz(:))';
orig=repmat(idx,size(onecube,1),1);

for i=1:size(onecube,2)
    face(:,i)=face(:,i)+orig(:);
end
face=unique(face,'rows');
face=mat2cell(face,ones(size(face,1),1));

if(nargout>=3)
    diffvec=cellfun(@diff,varargin,'UniformOutput',false);
    [xx,yy,zz]=ndgrid(diffvec{:});
    centroids=node(idx,:)+[xx(:) yy(:) zz(:)]*0.5;
end