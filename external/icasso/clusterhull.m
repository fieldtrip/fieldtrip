function [h_marker,txtcoord]=clusterhull(mode,coord,partition,color)
%function [h_marker,txtcoord]=clusterhull(mode,coord,partition,[color])
%
%PURPOSE
%
%To draw a cluster visualization where the objects (represented by
%points) belonging to the same cluster are presented by "cluster hulls",
%polygons where the edge of the polygon is the convex hull of the
%points that belong to the same cluster.
%
%INPUT
%
%%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
%mode      (string) 'edge' | 'fill'
%coord     (Mx2 matrix) coordinates of points (objects)
%partition (Mx1 vector) partition vector of objects (see
%            explanation for 'partition vector' in function hcluster)     
%[color]   (Kx3 matrix) each line is a 1x3 RGB vector that defines
%           color for each cluster. Default color is red for all 
%           clusters color(i,:)=[NaN NaN NaN] means that cluster(i)
%           is not drawn at all.
%
%OUTPUT
%
% h_marker (vector) graphic handles to all clusters. h_marker(i)=NaN if
%           color(i,:)=[NaN NaN NaN]
% txtcoord (Kx2 matrix) suggested coordinates for text labels.
%            txtcoord(i,:)=[NaN NaN] if cluster i is not drawn. 
%            Meaningful only in mode 'edge', returns always a
%            matrix of NaNs in mode 'fill'.
%
%DETAILS
%
%In mode 'edge' each cluster is represented by a convex hull if
%there are more than two points in the cluster. A two point cluster is
%represented by a line segment and a one point cluster by a small
%cross. The face of the convex hull is invisible.   
%
%In mode 'fill' the cluster hulls are filled with the specified
%color. Clusters having less than three members are ignored.
%
%NOTE The function always adds to a plot (turns 'hold on' temporarily).
%
%USED IN
%  icassoGraph

%COPYRIGHT NOTICE
%This function is a part of Icasso software library
%Copyright (C) 2003-2005 Johan Himberg
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% ver 1.2 100105 johan

holdStat=ishold;
hold on;

Ncluster=max(partition);

if nargin<4|isempty(color),
  color=repmat([1 0 0],Ncluster,1);
end

switch mode
 case 'edge'
  [h_marker,txtcoord]=clusteredge(coord,partition,color);
 case 'fill'
  h_marker=clusterfill(coord,partition,color);
  txtcoord=repmat(NaN,Ncluster,1);
 otherwise
  error('Unknown operation mode.');
end

% Some graphic settings

set(gca,'xtick',[],'ytick',[]); 
axis equal;

if ~holdStat,
  hold off;
end

function [h_marker,txtcoord]=clusteredge(coord,partition,color);

%% Number of clusters
Ncluster=max(partition);

%% Line width and +-sign size
LINEWIDTH=2;
MARKERSIZE=6;

for cluster=1:Ncluster,
  hullCoord=[];
  index=partition==cluster; Npoints=sum(index); 
  coord_=coord(index,:); 
  isVisible=all(isfinite(color(cluster,:)));
  if isVisible,
    if Npoints>=3,
      I=convhull(coord_(:,1),coord_(:,2)); 
      hullCoord=enlargehull([coord_(I,1) coord_(I,2)]);
      h_marker(cluster,1)=patch(hullCoord(:,1),hullCoord(:,2),0);
      set(h_marker(cluster,1),'edgecolor',color(cluster,:),'facecolor','none', ...
			'linewidth',LINEWIDTH);
    elseif size(coord_,1)==2,
      hullCoord(:,1)=[coord_(1,1) coord_(2,1)]'; hullCoord(:,2)=[coord_(1,2) coord_(2,2)]';
      h_marker(cluster,1)=plot(hullCoord(:,1)',hullCoord(:,2)','-');
      set(h_marker(cluster,1),'color',color(cluster,:),'linewidth',LINEWIDTH);
    else
      hullCoord(:,1)=coord_(1,1); hullCoord(:,2)=coord_(1,2);
      h_marker(cluster,1)=plot(coord_(1,1),coord_(1,2),'+');
      set(h_marker(cluster),'color',color(cluster,:),'markersize',MARKERSIZE);
  end
  txtcoord(cluster,:)=gettextcoord(hullCoord);
  else
    txtcoord(cluster,1:2)=NaN; h_marker(cluster,1)=NaN;
  end
end

function h_marker=clusterfill(coord,partition,color)


Ncluster=max(partition);

k=0;
for cluster=1:Ncluster,
  hullCoord=[];
  index=cluster==partition; Npoints=sum(index); 
  coord_=coord(index,:); 
  isVisible=all(isfinite(color(cluster,:)));
  if Npoints>=3 & isVisible,
    I=convhull(coord_(:,1),coord_(:,2)); 
    hullCoord=enlargehull([coord_(I,1) coord_(I,2)]);
    h_marker(cluster,1)=patch(hullCoord(:,1),hullCoord(:,2),0);
    set(h_marker(cluster),'edgecolor','none','facecolor',color(cluster,:));
  else
    h_marker(cluster,1)=NaN;
  end
end

function c=gettextcoord(coord)

[tmp,i]=min(coord(:,1));
c=coord(i,:);

function coord=enlargehull(coord)

S=1.1;

m=repmat(mean(coord(1:end-1,:),1),size(coord,1),1);

coord_=coord-m;

coord_=coord_.*S;
coord=coord_+m;

