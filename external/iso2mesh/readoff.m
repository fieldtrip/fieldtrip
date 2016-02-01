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
dim=sscanf(line,'OFF %d %d %d');
line=nonemptyline(fid);
if(size(dim,1)~=3)
    dim=sscanf(line,'%d',3);
    line=nonemptyline(fid);
end
nodalcount=3;
if(~isempty(line))
    [val nodalcount]=sscanf(line,'%f',inf);
else
    fclose(fid);
    return;
end
node=fscanf(fid,'%f',[nodalcount,dim(1)-1])';
node=[val(:)';node];

line=nonemptyline(fid);
facetcount=4;
if(~isempty(line))
    [val facetcount]=sscanf(line,'%f',inf);
else
    fclose(fid);
    return;
end
elem=fscanf(fid,'%f',[facetcount,dim(2)-1])';
elem=[val(:)';elem];
fclose(fid);
elem(:,1)=[];

if(size(elem,2)<=3)
    elem(:,1:3)=round(elem(:,1:3))+1;
else
    elem(:,1:4)=round(elem(:,1:4))+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str=nonemptyline(fid)
str='';
if(fid==0) error('invalid file'); end
while((isempty(regexp(str,'\S')) || ~isempty(regexp(str,'^#')))  && ~feof(fid))
    str=fgetl(fid);
    if(~ischar(str))
        str='';
        return;
    end
end
