function openedge=surfedge(f,varargin)
%
% openedge=surfedge(f)
%
% find the edge of an open surface or surface of a volume
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2007/11/21
%
% input:
%      f: input, surface face element list, dimension (be,3)
%
% output:
%      openedge: list of edges of the specified surface
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(isempty(f))
    openedge=[];
    return;
end

opt=varargin2struct(varargin{:});
findjunc=jsonopt('Junction',0,opt);

if(size(f,2)==3)
    edges=[f(:,[1,2]);
           f(:,[2,3]);
           f(:,[3,1])];             % create all the edges
elseif(size(f,2)==4)
    edges=[f(:,[1,2,3]);
           f(:,[2,1,4]);
           f(:,[1,3,4]);
           f(:,[2,4,3])];             % create all the edges
else
    error('surfedge only supports 2D and 3D elements');
end
% node4=[f(:,3);f(:,2);f(:,1)];   % node idx concatinated
edgesort=sort(edges,2);
[foo,ix,jx]=unique(edgesort,'rows');

if(isoctavemesh)
        u=unique(jx);
        if(size(f,2)==3 && findjunc)
            qx=u(hist(jx,u)>2);
        else
	    qx=u(hist(jx,u)==1);
        end
else
	vec=histc(jx,1:max(jx));
        if(size(f,2)==3 && findjunc)
            qx=find(vec>2);
        else
	    qx=find(vec==1);
        end
end
openedge=edges(ix(qx),:);
% node4=node4(ix(qx));
