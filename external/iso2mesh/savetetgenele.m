function savetetgenele(elem,fname)
%
% savetetgenele(elem,fname)
%
% save a mesh tetrahedral element list to tetgen .ele format
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%      elem: tetrahedral element list, dimension (ne,4)
%            columns beyound the 4rd column are treated as 
%            markers and attributes associated with the element
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
if(size(elem,2)>=6)
        hasprop=size(elem,2)-5;
        attrstr=repmat('%e ',1,hasprop);
end
if(size(elem,2)>=5)
        markers='%d';
end
elem(:,1:4)=elem(:,1:4)-1;
fprintf(fid,'%d %d %d\n',size(elem,1),4,hasprop+(size(elem,2)>=5));
fprintf(fid,['%d %d %d %d %d ' attrstr markers '\n'], [(1:size(elem,1))'-1 elem]');
fclose(fid);
