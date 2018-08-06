function [newnode,newelem,newface]=meshrefine(node,elem,varargin)
%
% [newnode,newelem,newface]=meshrefine(node,elem,face,opt)
%
% refine a tetrahedral mesh by adding new nodes or constraints
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input parameters:
%      node: existing tetrahedral mesh node list
%      elem: existing tetrahedral element list
%      face: (optional) existing tetrahedral mesh surface triangle list
%      opt:  options for mesh refinement:
%        if opt is a Nx3 array, opt is treated as a list of new nodes to
%          be inserted into the mesh. the new nodes must be located on the 
%          surface or inside the original mesh. external nodes are
%          discarded, unless the opt.extcmdopt is specified.
%        if opt is a vector with a length that equals to that of node,
%          it will be used to specify the desired edge-length at each node;
%          setting a node value to 0 will by-pass the refinement at this node
%        if opt is a vector with a length that equals to that of elem,
%          it will be used as the desired maximum element volume of each
%          tetrahedron; setting to 0 will by-pass the refinement of that element.
%        if opt is a struct, it can have the following fields:
%          opt.newnode: same as setting opt to an Nx3 array
%          opt.reratio: radius-edge ratio, by default, iso2mesh uses 1.414
%          opt.maxvol: maximum element volume
%          opt.sizefield: a vector specifying either the desired edge-length
%              at each node, or the maximum volume constraint within each 
%              tetrahedron, see above for details.
%          opt.extcmdopt: by default, meshrefine can only insert nodes
%              that are inside the original mesh. if one prefers to insert
%              nodes that are outside of the original mesh, one can define
%              this parameter to specify the meshing option (for tetgen)
%              for the extended domain, i.e. the convex hull including 
%              both the original and the external nodes. If not defined, 
%              '-Y' option is used by default (prevent tetgen from
%              inserting new nodes on the surface).
%          opt.extlabel: when external nodes are inserted, the new elements
%              will be assigned with an element label to group them
%              together, by default, this label is 0, unless opt.extlabel
%              is given
%          opt.extcorelabel: when external nodes are inserted, par of the 
%              new elements share the polyhedra between the inserted nodes,
%              these special elements will be marked by opt.extcorelabel, 
%              otherwise the label will be set to -1 
%
% outputs:
%      newnode: node coordinates of the tetrahedral mesh
%      newelem: element list of the tetrahedral mesh
%      newface: mesh surface element list of the tetrahedral mesh 
%             the last column denotes the boundary ID
%
% examples:
%
%     [node,face,elem]=meshasphere([0 0 0],24,1,2);
%     elem(:,5)=1;
%
%     % inserting nodes that are inside the original mesh
%     innernodes=double([1 1 1; 2 2 2; 3 3 3]);
%     [newno,newel]=meshrefine(node,elem,innernodes);
%     all(ismember(round(innernodes*1e10)*1e-10,round(newno*1e10)*1e-10,'rows'))
%     plotmesh(newno,[],newel,'x>-3')
%
%     % inserting nodes that are external to the original mesh
%     extnodes=double([-5 -5 25;-5 5 25;5 5 25;5 -5 25]);
%     [newno,newel]=meshrefine(node,elem,struct('newnode',extnodes,'extcmdopt','-Y'));
%     figure;
%     plotmesh(newno,[],newel,'x>-3')
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

exesuff=getexeext;
exesuff=fallbackexeext(exesuff,'tetgen');

newpt=[];
sizefield=[];
if(size(node,2)==4)
   sizefield=node(:,4);
   node=node(:,1:3);
end
opt=struct;
if(length(varargin)==1)
	face=[];
	if(isstruct(varargin{1}))
		opt=varargin{1};
        else
		if(length(varargin{1})==size(node,1) || length(varargin{1})==size(elem,1))
			sizefield=varargin{1};
		else
			newpt=varargin{1};
		end
	end
elseif(length(varargin)>=2)
    face=varargin{1};
    if(isstruct(varargin{2}))
        opt=varargin{2};
    else
        if(length(varargin{2})==size(node,1) || length(varargin{1})==size(elem,1))
                sizefield=varargin{2};
        else
                newpt=varargin{2};
        end
    end
else
	error('meshrefine requires at least 3 inputs');
end
if(isstruct(opt) && isfield(opt,'newnode'))
    newpt=opt.newnode;
end
if(isstruct(opt) && isfield(opt,'sizefield'))
    sizefield=opt.sizefield;
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

externalpt=[];
if(isstruct(opt) && isfield(opt,'extcmdopt'))
    isinside=tsearchn(node(:,1:3),elem(:,1:4),newpt(:,1:3));
    externalpt=newpt(isnan(isinside),:);
    newpt=newpt(find(isnan(isinside)==0),:);
end

if(~isempty(newpt))
	if(size(newpt,1)<4)
		newpt=[newpt; repmat(newpt(1,:),4-size(newpt,1),1)];
	end
	savetetgennode(newpt,mwpath('pre_refine.1.a.node'));
	moreopt=' -i ';
end

if(~isempty(sizefield))
    if(length(sizefield)==size(node,1))
        fid=fopen(mwpath('pre_refine.1.mtr'),'wt');
        fprintf(fid,'%d 1\n',size(sizefield,1));
        fprintf(fid,'%.16f\n',sizefield);
        fclose(fid);
	moreopt=[moreopt ' -qa '];
    else
        fid=fopen(mwpath('pre_refine.1.vol'),'wt');
        fprintf(fid,'%d\n',size(sizefield,1));
        fprintf(fid,'%d\t%.16f\n',[(1:size(sizefield,1))' sizefield]');
        fclose(fid);
	moreopt=[moreopt ' -qa '];
    end
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
elseif(~isempty(moreopt))
    system([' "' mcpath('tetgen') exesuff '" ' moreopt ' -r "' mwpath('pre_refine.1') '"']);
    [newnode,newelem,newface]=readtetgen(mwpath('pre_refine.2'));
else
    newnode=node;
    newelem=elem;
    newface=face;
end

if(~isempty(externalpt)) % user request to insert nodes that are outside of the original mesh
    % create a mesh including the external points
    externalpt=unique(externalpt,'rows');
    allnode=[newnode;externalpt];

    % define the convexhull as the external surface
    try
        outface=convhull(allnode,'simplify',false);
    catch
        outface=convhulln(allnode);
    end
    outface=sort(outface,2);
    face=volface(newelem(:,1:4));
    inface=sort(face(:,1:3),2);

    % define the surface that bounds the newly extended convex hull space
    bothsides=removedupelem([outface;inface]);

    % define a seed point to avoid meshing the interior space
    holelist=surfseeds(node,face(:,1:3));

    % mesh the extended space
    ISO2MESH_TETGENOPT=jsonopt('extcmdopt','-Y',opt);
    if(size(bothsides,1)>=size(inface,1))
        [no,el]=surf2mesh(allnode,bothsides,[],[],1,10,[],holelist);
    else
        [no,el]=surf2mesh(allnode,bothsides,[],[],1,10);
    end

    [isinside,map]=ismember(round(no*1e10)*1e-10,round(allnode*1e10)*1e-10,'rows');
    snid=[length(newnode)+1:length(allnode)];

    if(~isempty(map==0))
        oldsize=size(allnode,1);
        allnode=[allnode;no(map==0,:)];
        map(map==0)=oldsize+1:size(allnode,1);
    end
    % merge the external space with the original mesh
    el2=map(el(:,1:4));

    % label all new elements with -1
    if(size(newelem,2)==5)
        el2(:,5)=jsonopt('extlabel',0,opt);
        % search elements that contain source(s) and save their id(s)
        iselm=ismember(el2(:,1:4),snid);
        el2(sum(iselm,2)>=3,5)=jsonopt('extcorelabel',-1,opt);
    end

    % merge nodes/elements and replace the original ones
    newnode=allnode;
    newelem=[newelem;el2];
end

% read in the generated mesh

fprintf(1,'mesh refinement is complete\n');

