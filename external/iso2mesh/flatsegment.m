function mask=flatsegment(node,edge)
%
% mask=flatsegment(node,edge)
%
% decompose edge loops into flat segments alone arbitrary planes of the bounding box
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2008/04/08
%
% this code is fragile: it can not handle curves with many co-linear
% nodes near the corner point
%
% input:   
%    node:  x,y,z coordinates of each node of the mesh
%    edge:  input, a single vector separated by NaN, each segment
%           is a close-polygon consisted by node IDs 
%
% output:
%    mask:  output, a cell, each element is a close-polygon 
%           on x/y/z plane 
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

idx=edge;
nn=length(idx);
val=zeros(nn,1);
for i=1:nn
  tet=mod(i:i+3,nn);
  tet(find(tet==0))=nn;
  val(i)=(abs(det([node(idx(tet),:),ones(4,1)]))>1e-5);
end

val(end+1:end+2)=val(1:2);
mask={};
oldend=1;
for i=1:nn
	if(val(i)==1&val(i+1)==1&val(i+2)==0)
            val(i+2)=2;
            mask{count}=idx(oldend:i+2);
            count=count+1;
            oldend=i+2;
        else
            mask{count}=[idx(oldend:end);mask{1}];
            break;
        end
end
