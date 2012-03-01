function [node,elem,face]=readmedit(filename)
%
% [node,elem,face]=readmedit(filename)
%
% read Medit mesh format
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%    fname: name of the medit data file
%
% output:
%    node: node coordinates of the mesh
%    elem: list of elements of the mesh	    
%    face: list of surface triangles of the mesh	    
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

node=[];
elem=[];
face=[];
fid=fopen(filename,'rt');
while(~feof(fid))
    key=fscanf(fid,'%s',1);
    if(strcmp(key,'End')) break; end
    val=fscanf(fid,'%d',1);
    if(strcmp(key,'Vertices'))
        node=fscanf(fid,'%f',4*val);
        node=reshape(node,[4 val])';
    elseif(strcmp(key,'Triangles'))
        face=fscanf(fid,'%d',4*val);
        face=reshape(face,[4 val])';
    elseif(strcmp(key,'Tetrahedra'))
        elem=fscanf(fid,'%d',5*val);
        elem=reshape(elem,[5 val])';        
    end
end
fclose(fid);
