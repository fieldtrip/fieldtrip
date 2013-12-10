function [newnode,newelem,newface]=meshrefine(node,elem,varargin)
%
% [newnode,newelem,newface]=meshrefine(node,elem,face,opt)
%
% refine a tetrahedral mesh by adding new nodes or constraints
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input parameters:
%      node: existing tetrahedral mesh node list
%      elem: existing tetrahedral element list
%      face: (optional) existing tetrahedral mesh surface triangle list
%      opt:  options for mesh refinement:
%        if opt is a Nx3 array, opt is treated as a list of new nodes to
%          be inserted into the mesh (must be located on the surface or inside)
%        if opt is a struct, it can have the following fields:
%          opt.newnode: same as setting opt to an Nx3 array
%          opt.reratio: radius-edge ratio, by default, iso2mesh uses 1.414
%          opt.maxvol: maximum element volume
%
% outputs:
%      newnode: node coordinates of the tetrahedral mesh
%      newelem: element list of the tetrahedral mesh
%      newface: mesh surface element list of the tetrahedral mesh 
%             the last column denotes the boundary ID
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

exesuff=getexeext;
exesuff=fallbackexeext(exesuff,'tetgen');

newpt=[];
opt=struct;
if(length(varargin)==1)
	face=[];
	if(isstruct(varargin{1}))
		opt=varargin{1};
    else
		newpt=varargin{1};
	end
elseif(length(varargin)>=2)
    face=varargin{1};
    if(isstruct(varargin{2}))
        opt=varargin{2};
    else
        newpt=varargin{2};
    end
else
	error('meshrefine requires at least 3 inputs');
end
if(isstruct(opt) && isfield(opt,'newnode'))
    newpt=opt.newnode;
end

% call tetgen to create volumetric mesh
deletemeshfile(mwpath('pre_refine.*'));
deletemeshfile(mwpath('post_refine.*'));

moreopt='';
setquality=0;
if(isstruct(opt) && isfield(opt,'reratio'))
	moreopt=[moreopt sprintf(' -q %.10f ',opt.reratio)];
	setquality=1;
end
if(isstruct(opt) && isfield(opt,'maxvol'))
    moreopt=[moreopt sprintf(' -a %.10f ',opt.maxvol)];
end

if(~isempty(newpt))
	savetetgennode(newpt,mwpath('pre_refine.1.a.node'));
	moreopt=' -i ';
end
if(size(elem,2)==3 && setquality==0)
    if(~isempty(newpt))
        error('inserting new point can not be used for surfaces');
    end
    nedge=savegts(node, elem,mwpath('pre_refine.gts'));
    exesuff=fallbackexeext(getexeext,'gtsrefine');
elseif(size(elem,2)==3)
    savesurfpoly(node,elem,[],[],[],[],mwpath('pre_refine.poly'));
else
    savetetgennode(node, mwpath('pre_refine.1.node'));
    savetetgenele (elem, mwpath('pre_refine.1.ele'));
end

fprintf(1,'refining the input mesh ...\n');

if(size(elem,2)==3 && setquality==0)
    if(isstruct(opt) && isfield(opt,'scale'))
        moreopt=sprintf('%s -n %d ',moreopt,round(nedge*opt.scale));
    else
        error('you must give opt.scale value for refining a surface');
    end
end
if(isstruct(opt) && isfield(opt,'moreopt'))
	moreopt=[moreopt opt.moreopt];
end

if(size(elem,2)==3 && setquality==0)
    system([' "' mcpath('gtsrefine') exesuff '" ' moreopt ' < "' ...
          mwpath('pre_refine.gts') '" > "' mwpath('post_refine.gts') '"']);
    [newnode,newelem]=readgts(mwpath('post_refine.gts'));
    newface=newelem;
elseif(size(elem,2)==3)
    system([' "' mcpath('tetgen') exesuff '" ' moreopt ' -p -A "' mwpath('pre_refine.poly') '"']);
    [newnode,newelem,newface]=readtetgen(mwpath('pre_refine.1'));
else
    system([' "' mcpath('tetgen') exesuff '" ' moreopt ' -r "' mwpath('pre_refine.1') '"']);
    [newnode,newelem,newface]=readtetgen(mwpath('pre_refine.2'));
end

% read in the generated mesh

fprintf(1,'mesh refinement is complete\n');

