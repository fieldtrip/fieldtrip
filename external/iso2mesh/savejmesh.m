function savejmesh(node,face,elem,fname,varargin)
%
% savejmesh(node,face,elem,fname,opt)
%
% save a mesh to jMesh format
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2011/10/06
%
% input:
%      node: input, surface node list, dimension (nn,3)
%      face: input, surface face element list, dimension (be,3)
%      elem: input, tetrahedral element list, dimension (ne,4)
%      fname: output file name
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin==2)
   fname=face;
   face=[];
   elem=[];
end

if(nargin==3)
   fname=elem;
   elem=[];
end

fid=fopen(fname,'wt');
if(fid==-1)
    error('You do not have permission to save mesh files.');
end

mesh=struct('MeshNode',node,'MeshSurf',face,'MeshElem',elem,...
            'CreateTime',datestr(now),'Comment','Created by iso2mesh (http://iso2mesh.sf.net)');
fprintf(fid,'%s\n',savejson('',mesh,varargin{:}));

fclose(fid);
