function p=internalpoint(v,aloop)
%
% p=internalpoint(v,aloop)
%
% imperical function to find an internal point
% of a planar polygon
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2008/04/08
%
% input:   
%    v:     x,y,z coordinates of each node of the mesh
%    aloop:  input, a single vector separated by NaN, each segment
%             is a close-polygon consisted by node IDs 
% output:
%    p:   output, [x y z] of an internal point of aloop
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

p=[];
nd=v(aloop,:);
boxfacet=find(sum(abs(diff(nd)))<1e-2); % find a flat loop
if(length(boxfacet))   % if the loop is flat along x/y/z dir
     bf=boxfacet(1);    % no degeneracy allowed
     idx=setdiff([1 2 3],bf);
     
     p0=(nd(1,:)+nd(2,:))/2;
     pvec=complex(p0(idx(1)),p0(idx(2)));
     vec=nd(2,:)-nd(1,:);
     vec=complex(vec(idx(1)),vec(idx(2)))*exp(i*pi/2)*(1e-5)/sqrt(sum(vec.*vec));
     testpt=[real(pvec+vec) imag(pvec+vec);real(pvec-vec) imag(pvec-vec)];
     in=inpolygon(testpt(:,1),testpt(:,2), nd(:,idx(1)),nd(:,idx(2)));
     p=testpt(find(in>0),:);
     p([bf,idx(1),idx(2)])=[nd(1,bf),p];
end

if(length(p)==0|length(p)==2) 
    error('fail to find an internal point of curve');
end
