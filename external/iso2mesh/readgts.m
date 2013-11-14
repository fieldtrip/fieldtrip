function [node,elem,edges,edgemap]=readgts(fname)
%
% [node,elem,edges,edgemap]=readgts(fname)
%
% read GNU Triangulated Surface Format (GTS)
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2008/03/28
%
% input:
%    fname: name of the OFF data file
%
% output:
%    node: node coordinates of the mesh
%    elem: list of elements of the surface mesh
%    edges: the edge list section in the GTS file (optional)
%    edgemap: the face section (in terms of edge indices) in the GTS file
%             (optional)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

node=[];
elem=[];
fid=fopen(fname,'rt');
line=fgetl(fid);
dim=sscanf(line,'%d',3);
node   =fscanf(fid,'%f',[3,dim(1)])';
edges  =fscanf(fid,'%d',[2,dim(2)])';
edgemap=fscanf(fid,'%d',[3,dim(3)])';
fclose(fid);

edget=edges';
len=size(edgemap,1);
elem=reshape(edget(:,edgemap'),6,len)';
try
    for i=1:len
      elem(i,1:3)=unique(elem(i,:));
    end
catch
    error(sprint('invalid GTS face, id=%d\n',i));
end
elem=elem(:,1:3);