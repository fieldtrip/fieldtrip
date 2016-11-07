function [node,elem]=readsmf(fname)
%
% [node,elem]=readsmf(fname)
%
% read simple model format (SMF)
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2007/11/21
%
% input: 
%    fname: name of the	SMF data file
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
while(~feof(fid))
    line=fgetl(fid);
    if(line(1)=='v')
        dd=sscanf(line,'v %f %f %f');
        if(length(dd)==3)
            node=[node;dd];
        end
    elseif(line(1)=='f')
        dd=sscanf(line,'f %d %d %d');
        if(length(dd)==3)
            elem=[elem;dd];
        end
    end
end
fclose(fid);
node=reshape(node,3,length(node)/3)';
elem=reshape(elem,3,length(elem)/3)';
