function [coord_snapped] = warp_hermes2010(cfg, elec, surf)

% WARP_HERMES2010 projects the ECoG grid / strip onto a cortex hull
% using the algorithm described in Hermes et al. (2010,
% J Neurosci methods) in which electrodes are projected onto the pial
% surface using the orthogonal local norm vector to the grid. To align ECoG 
% electrodes to the pial surface, you first need to compute the cortex hull 
% with FT_PREPARE_MESH.
%
% See also FT_ELECTRODEREALIGN, FT_PREPARE_MESH, WARP_DYKSTRA2012

% Copyright (C) 2017, Dora Hermes, Arjen Stolk
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
%
% $Id$

disp('using warp algorithm described in Hermes et al. 2010: doi:10.1016/j.jneumeth.2009.10.005')

% determine grid dimensions (1st dim: number of arrays, 2nd dim: number of elecs in an array)
GridDim = determine_griddim(elec);

% grid/strip dimension-dependent settings:
% index -
% 0 fo global
% other for nr of electrodes for local
% checkdistance -
% 0, always project,
% 1 - if electrode within 2 mm of mask, do not project
% 2 - for strips, project all to closest distance
if isequal(GridDim(1),1) % strip
  index = 0; % 0 neighbors
  checkdistance = 2;
elseif isequal(GridDim(1),2) % two-lane strip
  index = 4; % 3 neighbors
  checkdistance = 1;
else % grid
  index = 5; % 4 neighbors
  checkdistance = 1;
end

checkdistance_dist=3; % 3 mm

if checkdistance==2
  disp('electrodes projected to closest point, no norm')
end

if mean(elec.elecpos(:,1))<0
  disp('left grid');
  % delete right hemisphere
  % surf.pos=surf.pos(surf.pos(:,1)<=0,:);
else
  disp('right grid');
  % delete right hemisphere
  % surf.pos=surf.pos(surf.pos(:,1)>=0,:);
end

if index==0 %global estimate of principal direction most orthogonal to array
  [v,d]=eig(cov(elec.elecpos)); %all vecs
  nm=v(:,find(diag(d)==min(diag(d)))); %vec we want
  nm=nm*sign(nm(1)*mean(elec.elecpos(:,1)));%check for left or rigth brain, invert nm if needed
end

out_ind=zeros(size(elec.elecpos,1),1);
for k=1:size(elec.elecpos,1)
  %sub array?
  if index>0, % get principal direction most orthogonal to sub-array
    [y,ind]=sort(Eucl_dist(elec.elecpos,elec.elecpos(k,:)'),'ascend');%select closest for sub-array
    [v,d]=eig(cov(elec.elecpos(ind(1:index),:))); %all vecs
    nm=v(:,find(diag(d)==min(diag(d)))); %vec we want
    nm=nm*sign(nm(1)*mean(elec.elecpos(:,1)));%check for left or right brain, invert nm if needed
  end
  %
  npls=[surf.pos(:,1)-elec.elecpos(k,1) surf.pos(:,2)-elec.elecpos(k,2) surf.pos(:,3)-elec.elecpos(k,3)]; %x,y,z lengths
  % calculate distance
  npls_dist=sqrt(sum((npls).^2,2));
  % check whether distance is < 3 mm
  distancesm2=0;
  if npls_dist(find(npls_dist==min(npls_dist)),:)<checkdistance_dist
    %disp(['distance < 3 mm electrode ' int2str(k) ]);
    distancesm2=1;
  end
  
  if checkdistance==1 && distancesm2==1 % electrode too close to surface to project
    out_ind(k)=find(npls_dist==min(npls_dist),1); %find minimum distance
  elseif checkdistance==2
    out_ind(k)=find(npls_dist==min(npls_dist),1); %find minimum distance
  else
    npls_unit=npls./repmat((sum(npls.^2,2).^.5),1,3); % normalize npls to get unit vector
    npdot=(npls_unit*nm); %distance along eigenvector direction (dot product)
    % only take surf.pos within 2.5 cm distance
    npdot(npls_dist>25)=0;
    %npdotrev=(npls_unit*-nm); % distance in reverse eigenvector direction
    out_ind(k)=find(abs(npdot)==max(abs(npdot)),1); %find minimum distance, max dot product
    %out_ind_rev(k)=find(npdotrev==max(npdotrev),1); %find minimum distance, max dot product
  end
end

coord_snapped=surf.pos(out_ind,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function distances = Eucl_dist(els, base)
distances = zeros(size(els,1),1);
for e = 1:size(els,1)
  distances(e) = sqrt(sum((els(e,:)-base').^2));
end
