function [cutpos,cutvalue,facedata,elemid]=qmeshcut(elem,node,value,cutat,varargin)
%
% [cutpos,cutvalue,facedata,elemid]=qmeshcut(elem,node,value,cutat)
%
% fast tetrahedral mesh slicer
%
% author:Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input: 
%   elem: integer array with dimensions of NE x 4, each row contains
%         the indices of all the nodes for each tetrahedron
%   node: node coordinates, 3 columns for x, y and z respectively
%   value: a scalar array with the length of node numbers, can have 
%          multiple columns 
%   cutat: cutat can have different forms:
%          if cutat is a 3x3 matrix, it defines a plane by 3 points: 
%                 cutat=[x1 y1 z1;x2 y2 z2;x3 y3 z3]
%          if cutat is a vector of 4 element, it defines a plane by
%                 a*x+b*y+c*z+d=0  and cutat=[a b c d]
%          if cutat is a single scalar, it defines an isosurface 
%                 inside the mesh at value=cutat
%          if cutat is a string, it defines an implicit surface
%                 at which the cut is made. it must has form expr1=expr2
%                 where expr1 expr2 are expressions made of x,y,z,v and
%                 constants
%          if cutat is a cell in the form of {'expression',scalar}, 
%                 the expression will be evaluated at each node to 
%                 produce a new value, then an isosurface is produced 
%                 at the levelset where new value=scalar; the 
%                 expression can contain constants and x,y,z,v
%
% output:
%   cutpos: all the intersections of mesh edges by the cutat
%           cutpos is similar to node, containing 3 columns for x/y/z
%   cutvalue: interpolated values at the intersections, with row number
%           being the num. of the intersections, column number being the 
%           same as "value".
%   facedata: define the intersection polygons in the form of patch "Faces"
%   elemid: the index of the elem in which each intersection polygon locates
%
%   without any output, qmeshcut generates a cross-section plot
%
% the outputs of this subroutine can be easily plotted using 
%
%  % if value has a length of node:
%     patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',cutvalue,'FaceColor','interp');
%
%  % if value has a length of elem:
%     patch('Vertices',cutpos,'Faces',facedata,'CData',cutvalue,'FaceColor','flat');
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

% get the coefficients of the cutat equation: ax+by+cz+d=0
if(nargin<4)
    error('qmeshcut requires at least 4 inputs');
end
if(size(value,1)~=size(node,1) && size(value,1)~=size(elem,1) && ~isempty(value))
    error('the length of value must be either that of node or elem');
end
if(isempty(value))
    cutvalue=[];
end
if(ischar(cutat) || (iscell(cutat) && length(cutat)==2 && ischar(cutat{1})))
    x=node(:,1);
    y=node(:,2);
    z=node(:,3);
    if(ischar(cutat))
        expr=regexp(cutat,'(.+)=(.+)','tokens','once'); %regexp(cutat,'=','split');
        if(length(expr)~=2) error('single expression must contain a single "=" sign'); end
        dist=eval(expr{1})-eval(expr{2});
    else
        dist=eval(cutat{1})-cutat{2};
    end
    if(all(dist<=0))
        asign=-double(dist<0);
        asign(asign==0)=1;
    else
        asign=double(dist>0);
        asign(asign==0)=-1;
    end
elseif(numel(cutat)==9 || numel(cutat)==4)
    if(numel(cutat)==9)
        [a,b,c,d]=getplanefrom3pt(cutat);
    else
        [a,b,c,d]=deal(cutat(:));
    end

    % compute which side of the cutat for all nodes in the mesh
    co=repmat([a b c],size(node,1),1);
    dist=sum( (co.*node)' )+d;
    asign=dist;
    asign(find(asign>=0))=1;
    asign(find(asign<0))=-1;
else
    if(size(value,1)~=size(node,1))
        error('must use nodal value list when cutting mesh at an isovalue');
    end
    dist=value-cutat;
    if(all(dist<=0))
        asign=-double(dist<0);
        asign(asign==0)=1;
    else
        asign=double(dist>0);
        asign(asign==0)=-1;
    end
end

% get all the edges of the mesh
esize=size(elem,2);
if(esize==4)
    edges=[elem(:,[1,2]);elem(:,[1,3]);elem(:,[1,4]);
           elem(:,[2,3]);elem(:,[2,4]);elem(:,[3,4])];
elseif(esize==3)
    edges=[elem(:,[1,2]);elem(:,[1,3]);elem(:,[2,3])];
elseif(esize==10)
    edges=[elem(:,[1,5]);elem(:,[1,8]);elem(:,[1,7]);
           elem(:,[2,5]);elem(:,[2,6]);elem(:,[2,9]);
	   elem(:,[3,6]);elem(:,[3,7]);elem(:,[3,10]);
	   elem(:,[4,8]);elem(:,[4,9]);elem(:,[4,10])];
end

% find all edges with two ends at the both sides of the plane
edgemask=sum(asign(edges),2);
cutedges=find(edgemask==0);
%edgemask=prod(asign(edges)');
%cutedges=find(edgemask<0);

% calculate the distances of the two nodes, and use them as interpolation weight 
cutweight=dist(edges(cutedges,:));
totalweight=diff(cutweight');

%caveat: if an edge is co-planar to the cutat, then totalweight will be 0
%        and dividing zero will cause trouble for cutweight

cutweight=abs(cutweight./repmat(totalweight(:),1,2));

% calculate the cross-cut position and the interpolated values

cutpos=node(edges(cutedges,1),:).*repmat(cutweight(:,2),[1 3])+...
       node(edges(cutedges,2),:).*repmat(cutweight(:,1),[1 3]);
if(size(value,1)==size(node,1))
  if(iscell(cutat) || ischar(cutat) || numel(cutat)==9 || numel(cutat)==4)
      cutvalue=value(edges(cutedges,1),:).*repmat(cutweight(:,2),[1 size(value,2)])+...
               value(edges(cutedges,2),:).*repmat(cutweight(:,1),[1 size(value,2)]);
  elseif(numel(cutat)==1)
      cutvalue=ones(size(cutpos,1),1)*cutat;
  end
end 
% organize all cross-cuts into patch facedata format

emap=zeros(size(edges,1),1);
emap(cutedges)=1:length(cutedges);
if(esize==10)
        emap=reshape(emap,[size(elem,1),12]); % 10-node element
else
	emap=reshape(emap,[size(elem,1),esize*(esize-1)/2]); % C^n_2
end
faceid=find(any(emap,2)==1);
facelen=length(faceid);

% cross-cuts can only be triangles or quadrilaterals for tetrahedral mesh
% (co-plannar mesh needs to be considered)

etag=sum(emap>0,2); % emap & etag are of length size(elem,1)

if(esize==3)  % surface mesh cut by a plane
	linecut=find(etag==2);
	lineseg=emap(linecut,:)';
	facedata=reshape(lineseg(find(lineseg)),[2,length(linecut)])';
    elemid=linecut;
    if(size(value,1)==size(elem,1) && ~exist('cutvalue','var'))
        cutvalue=value(elemid,:);
    end
	return;
end

tricut=find(etag==3);
quadcut=find(etag==4);
elemid=[tricut(:);quadcut(:)];
if(size(value,1)==size(elem,1) && ~exist('cutvalue','var'))
    cutvalue=value(elemid,:);
end
% fast way (vector-form) to get all triangles

tripatch=emap(tricut,:)';
tripatch=reshape(tripatch(find(tripatch)),[3,length(tricut)])';

% fast way to get all quadrilaterals in convexhull form (avoid using
% convhulln)

quadpatch=emap(quadcut,:)';
quadpatch=reshape(quadpatch(find(quadpatch)),[4,length(quadpatch)])';

% combine the two sets to create the final facedata
% using the matching-tetrahedra algorithm as shown in 
% https://visualization.hpc.mil/wiki/Marching_Tetrahedra

facedata=[tripatch(:,[1 2 3 3]); quadpatch(:,[1 2 4 3])];

% plot your results with the following command

if(nargout==0)
  patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',cutvalue,'facecolor','interp',varargin{:});
end
