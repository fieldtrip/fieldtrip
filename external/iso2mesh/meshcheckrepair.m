function [node,elem]=meshcheckrepair(node,elem,opt)
%
% [node,elem]=meshcheckrepair(node,elem,opt)
% 
% check and repair a surface mesh
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2008/10/10
%
% input/output:
%      node: input/output, surface node list, dimension (nn,3)
%      elem: input/output, surface face element list, dimension (be,3)
%      opt: options, including
%            'duplicated': remove duplicated elements
%            'isolated': remove isolated nodes
%            'deep': call external jmeshlib to remove non-manifold vertices
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin<3 || strcmp(opt,'duplicated'))
    l1=length(elem);
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
    if(length(eg)>0) 
        error('open surface found, you need to enclose it by padding zeros around the volume');
    end
end
exesuff=getexeext;
exesuff=fallbackexeext(exesuff,'meshfix');

if(nargin<3 || strcmp(opt,'deep'))
    deletemeshfile(mwpath('post_sclean.off'));
    saveoff(node,elem,mwpath('pre_sclean.off'));
    system([' "' mcpath('meshfix') exesuff '" "' mwpath('pre_sclean.off') '" "' mwpath('post_sclean.off') '"']);
    [node,elem]=readoff(mwpath('post_sclean.off'));
end
