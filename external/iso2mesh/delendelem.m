function elem=delendelem(elem,mask)
%
% elem=delendelem(elem,mask)
%
% delete elements whose nodes are all edge nodes
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2007/11/24
%
% input/output: 
%      elem: input/output, surface/volumetric element list
%      mask: of length of node number, =0 for internal nodes, =1 for edge nodes
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

badidx=sum(mask(elem)');
elem(find(badidx==size(elem,2)),:)=[];
