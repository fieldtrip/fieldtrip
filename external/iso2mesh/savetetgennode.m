function savetetgennode(node,fname)
%
% savetetgennode(node,fname)
%
% save a mesh node list to tetgen .node format
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%      node: node coordinates, dimension (nn,3)
%            columns beyound the 3rd column are treated as 
%            markers and attributes associated with the node
%      fname: output file name
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

hasprop=0;
attrstr='';
markers='';

fid=fopen(fname,'wt');
if(fid==0)
        error(['can not write to file ' fname]);
end
if(size(node,2)>=5)
        hasprop=size(node,2)-4;
        attrstr=repmat('%e ',1,hasprop);
end
if(size(node,2)>=4)
        markers='%d';
end
fprintf(fid,'%d %d %d %d\n',size(node,1),3,hasprop,size(node,2)>=4);
fprintf(fid,['%d %e %e %e ' attrstr markers '\n'], [(1:size(node,1))'-1 node]');
fclose(fid);
