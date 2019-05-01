function [node,elem]=readasc(fname)
%
% [node,elem]=readasc(fname)
%
% read FreeSurfer ASC mesh format
%
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
% date: 2009/04/02
% 
% input:
%      fname: name of the asc file
%
% output:
%      node: node positions of the mesh
%      elem: element list of the mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

node=[];
elem=[];
fid=fopen(fname,'rt');
if(fid==-1)
        error(['can not read file ' fname]);
end

line=fgetl(fid); % the first line is #!ascii ....
dim=fscanf(fid,'%d',2);
node=fscanf(fid,'%f',[4,dim(1)])';
elem=fscanf(fid,'%f',inf);
fclose(fid);

if(length(elem)==4*dim(2))
    elem=reshape(elem,[4,dim(2)])';
elseif(length(elem)==8*dim(2))
    elem=reshape(elem,[8,dim(2)])';
end

if(~any(node(:,end)))
	node=node(:,1:end-1);
end
if(~any(elem(:,end))) 
        elem=elem(:,1:end-1);
end

elem=elem+1;
