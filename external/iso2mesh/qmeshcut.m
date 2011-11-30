function [cutpos,cutvalue,facedata]=qmeshcut(elem,node,value,plane)
%
% [cutpos,cutvalue,facedata]=qmeshcut(elem,node,value,plane)
%
% fast tetrahedral mesh cross-section plot
%
% author:Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input: 
%   elem: integer array with dimensions of NE x 4, each row contains
%         the indices of all the nodes for each tetrahedron
%   node: node coordinates, 3 columns for x, y and z respectively
%   value: a scalar array with the length of node numbers, can have 
%          multiple columns 
%   plane: defines a plane by 3 points: plane=[x1 y1 z1;x2 y2 z2;x3 y3 z3]
%
% output:
%   cutpos: all the intersections of mesh edges by the plane
%           cutpos is similar to node, containing 3 columns for x/y/z
%   cutvalue: interpolated values at the intersections, with row number
%           being the num. of the intersections, column number being the 
%           same as "value".
%   facedata: define the intersection polygons in the form of patch "Faces"
%
% the outputs of this subroutine can be easily plotted using 
%     patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',cutvalue,...
%           'FaceColor','interp');
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

% get the coefficients of the plane equation: ax+by+cz+d=0
[a,b,c,d]=getplanefrom3pt(plane);

% compute which side of the plane for all nodes in the mesh
co=repmat([a b c],size(node,1),1);
dist=sum( (co.*node)' )+d;
asign=dist;
asign(find(asign>=0))=1;
asign(find(asign<0))=-1;

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

%caveat: if an edge is co-planar to the plane, then totalweight will be 0
%        and dividing zero will cause trouble for cutweight

cutweight=abs(cutweight./repmat(totalweight(:),1,2));

% calculate the cross-cut position and the interpolated values

cutpos=node(edges(cutedges,1),:).*repmat(cutweight(:,2),[1 3])+...
       node(edges(cutedges,2),:).*repmat(cutweight(:,1),[1 3]);
cutvalue=value(edges(cutedges,1),:).*repmat(cutweight(:,2),[1 size(value,2)])+...
       value(edges(cutedges,2),:).*repmat(cutweight(:,1),[1 size(value,2)]);
   
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

etag=sum(emap>0,2);

if(esize==3)  % surface mesh cut by a plane
	linecut=find(etag==2);
	lineseg=emap(linecut,:)';
	facedata=reshape(lineseg(find(lineseg)),[2,length(linecut)])';
	return;
end

tricut=find(etag==3);
quadcut=find(etag==4);

% fast way (vector-form) to get all triangles

tripatch=emap(tricut,:)';
tripatch=reshape(tripatch(find(tripatch)),[3,length(tricut)])';

% fast wall to get all quadrilaterals in convexhull form (avoid using convhulln)

quadpatch=emap(quadcut,:)';
quadpatch=reshape(quadpatch(find(quadpatch)),[4,length(quadpatch)])';

% combine the two sets to create the final facedata
% using the matching-tetrahedra algorithm as shown in 
% https://visualization.hpc.mil/wiki/Marching_Tetrahedra

facedata=[tripatch(:,[1 2 3 3]); quadpatch(:,[1 2 4 3])];

% plot your results with the following command

%patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',cutvalue,'facecolor','interp');

