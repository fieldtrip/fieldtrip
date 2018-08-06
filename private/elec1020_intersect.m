function [cnt1, cnt2] = elec1020_intersect(pnt, dhk, v1, v2, v3, flag);

% ELEC1020_INTERSECT determines the intersection of a mesh with a plane

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

% Copyright (C) 2002 Robert Oostenveld

if nargin<6
  flag = 0;
end

npnt = size(pnt, 1);
ndhk = size(dhk, 1);

for i=1:npnt
  % determine on which side of the plane each vertex lies
  pside(i) = ptriside(v1, v2, v3, pnt(i,:));
end
    
for i=1:ndhk
  % find the triangles that are intersected by the plane
  dcut(i) = abs(sum(pside(dhk(i,:))))<2;
end 

% further computations will only be done on those selected triangles
sel = dhk(find(dcut),:);
nsel = length(find(dcut));

cnt1 = [];
cnt2 = [];
for i=1:nsel
  [l1, l2] = tritrisect(v1, v2, v3, pnt(sel(i,1),:), pnt(sel(i,2),:), pnt(sel(i,3),:));
  cnt1 = [cnt1; l1];
  cnt2 = [cnt2; l2];
end

% determine which segments connect to each other
indx = 1;
for i=1:nsel
  for j=setdiff(1:nsel, indx)
    d1 = pntdist(cnt2(i,:), cnt1(j,:));
    d2 = pntdist(cnt2(i,:), cnt2(j,:));
    if d1<eps || d2<eps
      indx = [indx; j]; 
    end
  end
end

cnt1 = cnt1(indx, :);
cnt2 = cnt2(indx, :);

% flip segments if they are the wrong way around
for i=1:(size(cnt1,1)-1)
  d1 = pntdist(cnt2(i,:), cnt1(i+1,:));
  d2 = pntdist(cnt2(i,:), cnt2(i+1,:));
  if d2<d1
    tmp1 = cnt1(i+1,:);
    tmp2 = cnt2(i+1,:);
    cnt1(i+1,:) = tmp2;
    cnt2(i+1,:) = tmp1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag==1

  e1 = (v3-v1)/norm(v3-v1);
  e2 = (v2-v1)/norm(v2-v1);
  e3 = cross(e2, e1);

  t1 = v1;
  t2 = v1 + e1;
  t3 = v1 + e3;

  ncnt=size(cnt1,1);

  new_cnt1 = [];
  new_cnt2 = [];

  for i=1:ncnt
    side(1) = ptriside(t1, t2, t3, cnt1(i,:)); 
    side(2) = ptriside(t1, t2, t3, cnt2(i,:)); 
    if all(side==1)
      new_cnt1 = [new_cnt1; cnt1(i,:)];
      new_cnt2 = [new_cnt2; cnt2(i,:)];
    elseif any(side==0) && any(side==1)
      new_cnt1 = [new_cnt1; cnt1(i,:)];
      new_cnt2 = [new_cnt2; cnt2(i,:)];
    elseif any(side==1) && any(side==-1)
      if side==-1
        new_cnt1 = [new_cnt1; ptriproj(t1, t2, t3, cnt1(i,:))];
        new_cnt2 = [new_cnt2; cnt2(i,:)];
      else
        new_cnt1 = [new_cnt1; cnt1(i,:)];
        new_cnt2 = [new_cnt2; ptriproj(t1, t2, t3, cnt2(i,:))];
      end
    end
  end

  cnt1 = new_cnt1;
  cnt2 = new_cnt2;

end

