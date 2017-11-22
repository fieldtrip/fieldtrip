function hh=plotedges(node,edges,varargin)
%
% hm=plotedges(node,edges,opt)
%   or
% hm=plotedges(node,loops,opt)
%
% plot a 3D polyline or close loop (1d manifold)
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input: 
%      node: node coordinates, dimension (nn,3); if node has a 
%            4th column, it will be used to set the color at each node.
%      edges:edge list: a 2-column index array, with each row being
%            an edge connecting the two indexed node
%      loops:loops is an NaN separated integer array, with each segment
%            denoting a 3D polyline or loop represented by a list of node
%            indices
%      opt:  additional options for the plotting, see plotmesh
%
% output:
%      hm: handle or handles (vector) to the plotted surfaces
%
% example:
%
%   h=plotedges(node,[1 2 3 4 5 nan 6 7 8 9]);
%   h=plotedges(node,edges,'marker','o','linewidth',2,'color','r');
% 
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

edlen=size(edges,1);
hh=[];
if(isempty(edges))
    return;
end

rngstate = rand ('state');

if(size(edges,1)==1 || size(edges,2)==1) % a loop: single column/row edges
    randseed=hex2dec('623F9A9E'); % "U+623F U+9A9E"
    if(isoctavemesh) randseed=randseed+3; end
    if(~isempty(getvarfrom({'caller','base'},'ISO2MESH_RANDSEED')))
        randseed=getvarfrom({'caller','base'},'ISO2MESH_RANDSEED');
    end
    rand('state',randseed);

    loops=edges(:)';
    if(~isnan(loops(end)))
        loops(end+1)=nan;
    end
    seg=find(isnan(loops));
    seglen=length(seg);
    seghead=1;
    for i=1:seglen
        h=plotmesh(node(loops(seghead:seg(i)-1),:), 'color',rand(3,1), varargin{:});
        hh=[hh h];
        seghead=seg(i)+1;
    end
else  % an edge list
 if(size(node,2)==2) % 2D polyline
    for i=1:edlen
       h = line([node(edges(i,1),1) node(edges(i,2),1)],[node(edges(i,1),2) node(edges(i,2),2)], varargin{:});
       hh=[hh,h];
    end
 else  % 3D polyline
    for i=1:edlen
       h = line([node(edges(i,1),1) node(edges(i,2),1)],[node(edges(i,1),2) ...
              node(edges(i,2),2)],[node(edges(i,1),3) node(edges(i,2),3)], varargin{:});
       hh=[hh,h];
    end
 end
end

rand ('state',rngstate);
