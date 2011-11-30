function [node,elem]=readoff(fname)
%
% [node,elem]=readoff(fname)
%
% read Geomview Object File Format (OFF)
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2008/03/28
%
% input:
%    fname: name of the OFF data file
%
% output:
%    node: node coordinates of the mesh
%    elem: list of elements of the mesh	    
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

node=[];
elem=[];
fid=fopen(fname,'rt');
line=fgetl(fid);
dim=fscanf(fid,'%d',3);
node=fscanf(fid,'%f',[3,dim(1)])';
elem=fscanf(fid,'%f',inf);
fclose(fid);
if(length(elem)==4*dim(2))
    elem=reshape(elem,[4,dim(2)])';
elseif(length(elem)==8*dim(2))
    elem=reshape(elem,[8,dim(2)])';
end
if(size(elem,2)<=3)
    elem=round(elem(:,2:3))+1;
else
    elem=round(elem(:,2:4))+1;
end
