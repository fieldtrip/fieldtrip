function [node,elem,face]=readtetgen(fstub)
%
% [node,elem,face]=readtetgen(fstub)
%
% read tetgen output files
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2007/11/21
%
% input:
%    fstub: file name stub
%
% output:
%    node: node coordinates of the tetgen mesh
%    elem: tetrahedra element list of the tetgen mesh
%    face: surface triangles of the tetgen mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

% read node file
fp=fopen([fstub,'.node'],'rt');
if(fp==0) 
	error('node file is missing!'); 
end
[dim,count] = fscanf(fp,'%d',4);
if(count<4) error('wrong node file'); end
node=fscanf(fp,'%f',[4,dim(1)]);
idx=node(1,:);
node=node(2:4,:)';
fclose(fp);

% read element file
fp=fopen([fstub,'.ele'],'rt');
if(fp==0) 
        error('elem file is missing!'); 
end
[dim,count] = fscanf(fp,'%d',3);
if(count<3) error('wrong elem file'); end
elem=fscanf(fp,'%d',[dim(2)+dim(3)+1,dim(1)]);
elem=elem';
elem(:,1)=[];
elem(:,1:dim(2))=elem(:,1:dim(2))+(1-idx(1));
fclose(fp);

% read surface mesh file
fp=fopen([fstub,'.face'],'rt');
if(fp==0)
        error('surface data file is missing!');
end
[dim,count] = fscanf(fp,'%d',2);
if(count<2) error('wrong surface file'); end
face=fscanf(fp,'%d',[5,dim(1)]);
face=[face(2:end-1,:)+1;face(end,:)]';
fclose(fp);

