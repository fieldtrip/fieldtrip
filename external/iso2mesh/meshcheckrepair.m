function [node,elem]=meshcheckrepair(node,elem,opt,varargin)
%
% [node,elem]=meshcheckrepair(node,elem,opt)
% 
% check and repair a surface mesh
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2008/10/10
%
% input/output:
%      node: input/output, surface node list, dimension (nn,3)
%      elem: input/output, surface face element list, dimension (be,3)
%      opt: options, including
%            'dupnode': remove duplicated nodes
%            'dupelem' or 'duplicated': remove duplicated elements
%            'dup': both above
%            'isolated': remove isolated nodes
%            'open': abort when open surface is found
%            'deep': call external jmeshlib to remove non-manifold vertices
%            'meshfix': repair a closed surface by the meshfix utility (new)
%                       it can remove self-intersecting elements and fill holes
%            'intersect': test a surface for self-intersecting elements
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin<3 || strcmp(opt,'dupnode')|| strcmp(opt,'dup'))
    l1=size(node,1);
    [node,elem]=removedupnodes(node,elem);
    l2=size(node,1);
    if(l2~=l1) fprintf(1,'%d duplicated nodes were removed\n',l1-l2); end
end

if(nargin<3 || strcmp(opt,'duplicated')|| strcmp(opt,'dupelem')|| strcmp(opt,'dup'))
    l1=size(elem,1);
    elem=removedupelem(elem);
    l2=length(elem);
    if(l2~=l1) fprintf(1,'%d duplicated elements were removed\n',l1-l2); end
end

if(nargin<3 || strcmp(opt,'isolated'))
    l1=length(node);
    [node,elem]=removeisolatednode(node,elem);
    l2=length(node);
    if(l2~=l1) fprintf(1,'%d isolated nodes were removed\n',l1-l2); end
end

if(nargin==3 && strcmp(opt,'open'))
    eg=surfedge(elem);
    if(~isempty(eg)) 
        error('open surface found, you need to enclose it by padding zeros around the volume');
    end
end

if(nargin<3 || strcmp(opt,'deep'))
    exesuff=getexeext;
    exesuff=fallbackexeext(exesuff,'jmeshlib');
    deletemeshfile(mwpath('post_sclean.off'));
    saveoff(node(:,1:3),elem(:,1:3),mwpath('pre_sclean.off'));
    system([' "' mcpath('jmeshlib') exesuff '" "' mwpath('pre_sclean.off') '" "' mwpath('post_sclean.off') '"']);
    [node,elem]=readoff(mwpath('post_sclean.off'));
end

exesuff=fallbackexeext(getexeext,'meshfix');
extra=varargin2struct(varargin{:});
moreopt=' -q -a 0.01 ';
if(isstruct(extra) && isfield(extra,'MeshfixParam'))
    moreopt=extra.MeshfixParam;
end

if(nargin>=3 && strcmp(opt,'meshfix'))
    deletemeshfile(mwpath('pre_sclean.off'));
    deletemeshfile(mwpath('pre_sclean_fixed.off'));
    saveoff(node,elem,mwpath('pre_sclean.off'));
    system([' "' mcpath('meshfix') exesuff '" "' mwpath('pre_sclean.off') ...
        '" ' moreopt]);
    [node,elem]=readoff(mwpath('pre_sclean_fixed.off'));
end

if(nargin>=3 && strcmp(opt,'intersect'))
    moreopt=sprintf(' -q --no-clean --intersect -o "%s"',mwpath('pre_sclean_inter.msh'));
    deletemeshfile(mwpath('pre_sclean.off'));
    deletemeshfile(mwpath('pre_sclean_inter.msh'));
    saveoff(node,elem,mwpath('pre_sclean.off'));
    system([' "' mcpath('meshfix') exesuff '" "' mwpath('pre_sclean.off') ...
        '" ' moreopt]);
    %[node,elem]=readoff(mwpath('pre_sclean_inter.off'));
end
