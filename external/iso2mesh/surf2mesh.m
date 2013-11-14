function [node,elem,face]=surf2mesh(v,f,p0,p1,keepratio,maxvol,regions,holes,forcebox)
%
% [node,elem,face]=surf2mesh(v,f,p0,p1,keepratio,maxvol,regions,holes,forcebox)
%
% create quality volumetric mesh from isosurface patches
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/11/24
%
% input parameters:
%      v: input, isosurface node list, dimension (nn,3)
%         if v has 4 columns, the last column specifies mesh density near each node
%      f: input, isosurface face element list, dimension (be,3)
%      p0: input, coordinates of one corner of the bounding box, p0=[x0 y0 z0]
%      p1: input, coordinates of the other corner of the bounding box, p1=[x1 y1 z1]
%      keepratio: input, percentage of elements being kept after the simplification
%      maxvol: input, maximum tetrahedra element volume
%      regions: list of regions, specifying by an internal point for each region
%      holes: list of holes, similar to regions
%      forcebox: 1: add bounding box, 0: automatic
%
% outputs:
%      node: output, node coordinates of the tetrahedral mesh
%      elem: output, element list of the tetrahedral mesh
%      face: output, mesh surface element list of the tetrahedral mesh 
%             the last column denotes the boundary ID
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fprintf(1,'generating tetrahedral mesh from closed surfaces ...\n');

exesuff=getexeext;
exesuff=fallbackexeext(exesuff,'tetgen');

if(keepratio>1 | keepratio<0)
   warn(['The "keepratio" parameter is required to be between 0 and 1. '...
         'Your input is out of this range. surf2mesh will not perform '...
	 'simplification. Please double check to correct this.']);
end

% first, resample the surface mesh with cgal
if(keepratio<1-1e-9 & ~iscell(f))
	fprintf(1,'resampling surface mesh ...\n');
	[no,el]=meshresample(v(:,1:3),f(:,1:3),keepratio);
	el=unique(sort(el,2),'rows');

	% then smooth the resampled surface mesh (Laplace smoothing)

	%% edges=surfedge(el);  % disable on 12/05/08, very slow on octave
	%% mask=zeros(size(no,1),1);
	%% mask(unique(edges(:)))=1;  % =1 for edge nodes, =0 otherwise
	%[conn,connnum,count]=meshconn(el,length(no));
	%no=smoothsurf(no,mask,conn,2);

	% remove end elements (all nodes are edge nodes)
	%el=delendelem(el,mask);
else
	no=v;
	el=f;
end
if(nargin==6)
	regions=[];
	holes=[];
elseif(nargin==7)
	holes=[];
end

dobbx=0;
if(nargin>=9)
	dobbx=forcebox;
end

% dump surface mesh to .poly file format
if(~iscell(el) & ~isempty(no) & ~isempty(el))
	saveoff(no(:,1:3),el(:,1:3),mwpath('post_vmesh.off'));
end
deletemeshfile(mwpath('post_vmesh.mtr'));
savesurfpoly(no,el,holes,regions,p0,p1,mwpath('post_vmesh.poly'),dobbx);

moreopt='';
if(size(no,2)==4)
   moreopt=[moreopt ' -m '];
end
% call tetgen to create volumetric mesh
deletemeshfile(mwpath('post_vmesh.1.*'));
fprintf(1,'creating volumetric mesh from a surface mesh ...\n');

fprintf(1,sprintf('\n%s\n%s\n',...
     'WARNING: the license for "tetgen" is non-free and does not permit commercial use.', ...
     'Please use the "cgalmesh" or "cgalpoly" options where free-software is desired.'));

cmdopt=getvarfrom({'caller','base'},'ISO2MESH_TETGENOPT');
if(isempty(cmdopt))
  system([' "' mcpath('tetgen') exesuff '" -A -q1.414a' num2str(maxvol) ' ' moreopt ' "' mwpath('post_vmesh.poly') '"']);
else
  system([' "' mcpath('tetgen') exesuff '" ' cmdopt ' "' mwpath('post_vmesh.poly') '"']);
end

% read in the generated mesh
[node,elem,face]=readtetgen(mwpath('post_vmesh.1'));

fprintf(1,'volume mesh generation is complete\n');

