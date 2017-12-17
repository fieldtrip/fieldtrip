function [cnt1, cnt2] = elec1020_follow(pnt, dhk, v1, v2, v3, feedback)

% ELEC1020_FOLLOW

% Copyright (C) 2003, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.

if nargin<6
  feedback = false;
end

% compute the distribution of edge lengths
edge = [
  pnt(dhk(:,1),:) - pnt(dhk(:,2),:)
  pnt(dhk(:,2),:) - pnt(dhk(:,3),:)
  pnt(dhk(:,3),:) - pnt(dhk(:,1),:)
  ];
edgelength = sqrt(sum(edge.^2,2));

% the tolerance depends on the edge length
tolerance       = min(edgelength)*1e-6;
tolerance_limit = min(edgelength)*1e-3;

npnt = size(pnt,1);
ndhk = size(dhk,1);

pside = nan(1,npnt);
for i=1:npnt
  % determine on which side of the plane each vertex lies
  pside(i) = ptriside(v1, v2, v3, pnt(i,:), tolerance);
end

% dcut = zeros(ndhk,1);
% for i=1:ndhk
%   % find the triangles that are intersected by the plane
%   if sum(pside(dhk(i,:))==0)==2
%     dcut(i) = 1;
%   elseif sum(pside(dhk(i,:))==0)==1 & sum(pside(dhk(i,:))==1)==1 & sum(pside(dhk(i,:))==-1)==1
%     dcut(i) = 1;
%   elseif sum(pside(dhk(i,:))==1)==2 & sum(pside(dhk(i,:))==-1)==1
%     dcut(i) = 1;
%   elseif sum(pside(dhk(i,:))==1)==1 & sum(pside(dhk(i,:))==-1)==2
%     dcut(i) = 1;
%   end
% end
tmp  = pside(dhk);
dcut = true(ndhk,1);
dcut(all(tmp== 1,2)) = false;
dcut(all(tmp==-1,2)) = false;

% continue working with only the intersecting triangles
dhk = dhk(dcut,:);
ndhk = size(dhk,1);

% for each triangle determine the neighbouring triangles
neighb = zeros(ndhk,ndhk);
for i=1:ndhk
  for j=(i+1):ndhk
    if dhk(i,1)==dhk(j,1)
      neighb(i,j) = 1;
    elseif dhk(i,1)==dhk(j,2)
      neighb(i,j) = 1;
    elseif dhk(i,1)==dhk(j,3)
      neighb(i,j) = 1;
    elseif dhk(i,2)==dhk(j,1)
      neighb(i,j) = 1;
    elseif dhk(i,2)==dhk(j,2)
      neighb(i,j) = 1;
    elseif dhk(i,2)==dhk(j,3)
      neighb(i,j) = 1;
    elseif dhk(i,3)==dhk(j,1)
      neighb(i,j) = 1;
    elseif dhk(i,3)==dhk(j,2)
      neighb(i,j) = 1;
    elseif dhk(i,3)==dhk(j,3)
      neighb(i,j) = 1;
    else
      neighb(i,j) = 0;
    end
    neighb(j,i) = neighb(i,j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a contour from v1 via v2 to v3

% find the nearest triangle on which point v1 projects
v1_dist = inf;
for i=1:ndhk
  % shift a fraction towards v2 to avoid starting in the vertex of a wrong triangle
  [proj, dist] = ptriproj(pnt(dhk(i,1),:), pnt(dhk(i,2),:), pnt(dhk(i,3),:), v1+tolerance*(v2-v1)/norm(v2-v1), 1);
  if dist<v1_dist
    v1_dist = dist;
    v1_indx = i;
    v1_proj = ptriproj(pnt(dhk(i,1),:), pnt(dhk(i,2),:), pnt(dhk(i,3),:), v1, 1);
  end
end

% find the nearest triangle on which point v3 projects
v3_dist = inf;
for i=1:ndhk
  [proj, dist] = ptriproj(pnt(dhk(i,1),:), pnt(dhk(i,2),:), pnt(dhk(i,3),:), v3, 1);
  if dist<v3_dist
    v3_dist = dist;
    v3_indx = i;
    v3_proj = proj;
  end
end

% intersect the triangle containing the projection of v1 with the plane
[l1, l2] = tritrisect(v1, v2, v3, pnt(dhk(v1_indx,1),:), pnt(dhk(v1_indx,2),:), pnt(dhk(v1_indx,3),:));
if pntdist(l1,v2) < pntdist(l2,v2)
  % remember the projection point and l1 as the first line segment
  cnt1 = [v1_proj];
  cnt2 = [l1];
  prev_proj = l1;
  prev_indx = v1_indx;
else
  % remember the projection point and l2 as the first line segment
  cnt1 = [v1_proj];
  cnt2 = [l2];
  prev_proj = l2;
  prev_indx = v1_indx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while(1)
  
  for i=find(neighb(prev_indx,:))
    
    ncnt = size(cnt1,1);
    c1 = pnt(dhk(i,1),:);
    c2 = pnt(dhk(i,2),:);
    c3 = pnt(dhk(i,3),:);
%     
%     if feedback
%       ft_plot_mesh(struct('pos', pnt, 'tri', dhk(i,:)), 'edgecolor', 'b', 'facecolor', 'none')
%       ft_plot_mesh(c1, 'vertexsize', 30, 'vertexcolor', 'r')
%       ft_plot_mesh(c2, 'vertexsize', 30, 'vertexcolor', 'r')
%       ft_plot_mesh(c3, 'vertexsize', 30, 'vertexcolor', 'r')
%       ft_plot_mesh(prev_proj, 'vertexsize', 30, 'vertexcolor', 'g')
%       drawnow
%     end
    
    [proj, dist] = ptriproj(c1, c2, c3, prev_proj, 1);
    
    if dist<tolerance
      [l1, l2] = tritrisect(v1, v2, v3, c1, c2, c3);
      
      if pntdist(l1, cnt1(ncnt,:))<tolerance && pntdist(l2, cnt2(ncnt,:))<tolerance
        continue
      elseif pntdist(l1, cnt2(ncnt,:))<tolerance && pntdist(l2, cnt1(ncnt,:))<tolerance
        continue
      end
      
      if pntdist(l1, prev_proj) < pntdist(l2, prev_proj)
        cnt1 = [cnt1; l1];
        cnt2 = [cnt2; l2];
        prev_proj = l2;
        prev_indx = i;
        break
      else
        cnt1 = [cnt1; l2];
        cnt2 = [cnt2; l1];
        prev_proj = l1;
        prev_indx = i;
        break
      end
    end
    
  end
  
  % stop if no new segment was added
  if ~(size(cnt1,1)>ncnt)
    tolerance = 2*tolerance;
    if tolerance>=tolerance_limit
      ft_warning('premature end of contour')
      break
    else
      ft_warning('increasing tolerance');
    end
  end
  
  % stop if we arrive on the triangle with the endpoint
  if pntdist(prev_proj, v3_proj)<tolerance
    % replace the current endpoint with the projection of v3
    cnt2(size(cnt2,1),:) = v3_proj;
    break
  end
  
  % stop if we arrive on the triangle with the endpoint
  if prev_indx==v3_indx
    % replace the current endpoint with the projection of v3
    cnt2(size(cnt2,1),:) = v3_proj;
    break
  end
  
end % while

if feedback
  X = [cnt1(:,1) cnt2(:,1)];
  Y = [cnt1(:,2) cnt2(:,2)];
  Z = [cnt1(:,3) cnt2(:,3)];
  line(X, Y, Z, 'color', 'k')
  drawnow
end

end % function
