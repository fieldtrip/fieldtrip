function [coord_snapped] = warp_hermes2010(cfg, elec, surf)

% WARP_HERMES2012 projects the ECoG grid / strip onto a cortex hull
% while minimizing the distance from original positions and the
% deformation of the grid. To align ECoG electrodes to the pial surface,
% you first need to compute the cortex hull with FT_PREPARE_MESH.
%
% WARP_HERMES2012 uses the algorithm described in Hermes et al. (2010,
% J Neurosci methods) in which electrodes are projected onto the pial 
% surface using the orthogonal local norm vector to the grid
%
% See also FT_ELECTRODEREALIGN, FT_PREPARE_MESH

% Copyright (C) 2010-2017, Dora Hermes, Arjen Stolk
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

els = elec.elecpos;
gs = surf.pos;

% grid/strip dimension-dependent settings:
% index - 
    % 0 fo global
    % other for nr of electrodes for local 
% checkdistance - 
    % 0, always project, 
    % 1 - if electrode within 2 mm of mask, do not project
    % 2 - for strips, project all to closest distance
GridDim = determine_griddim(elec);
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
    
lcl_num=index;
checkdistance_dist=3; % 3 mm

if checkdistance==2
    disp('electrodes projected to closest point, no norm')
end
    
if mean(els(:,1))<0
    disp('left grid');
    % delete right hemisphere
    % gs=gs(gs(:,1)<=0,:);
else
    disp('right grid');
    % delete right hemisphere
    % gs=gs(gs(:,1)>=0,:);
end

if lcl_num==0 %global estimate of principal direction most orthogonal to array
    [v,d]=eig(cov(els)); %all vecs
    nm=v(:,find(diag(d)==min(diag(d)))); %vec we want
    nm=nm*sign(nm(1)*mean(els(:,1)));%check for left or rigth brain, invert nm if needed
end

out_ind=zeros(size(els,1),1);
out_ind_rev=zeros(size(els,1),1);
for k=1:size(els,1)
    %sub array?
    if lcl_num>0, % get principal direction most orthogonal to sub-array
        [y,ind]=sort(dist_hermes(els,els(k,:)'),'ascend');%select closest for sub-array
        [v,d]=eig(cov(els(ind(1:lcl_num),:))); %all vecs
        nm=v(:,find(diag(d)==min(diag(d)))); %vec we want
        nm=nm*sign(nm(1)*mean(els(:,1)));%check for left or right brain, invert nm if needed
    end
    %
    npls=[gs(:,1)-els(k,1) gs(:,2)-els(k,2) gs(:,3)-els(k,3)]; %x,y,z lengths
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
        % only take gs within 2.5 cm distance
        npdot(npls_dist>25)=0;
        %npdotrev=(npls_unit*-nm); % distance in reverse eigenvector direction
        out_ind(k)=find(abs(npdot)==max(abs(npdot)),1); %find minimum distance, max dot product
        %out_ind_rev(k)=find(npdotrev==max(npdotrev),1); %find minimum distance, max dot product
    end
end

coord_snapped=gs(out_ind,:);
% out_els_rev=gs(out_ind_rev,:);
% % check distance to new electrodes
% out_els_dist=sqrt(sum((els-out_els).^2,2));
% out_els_rev_dist=sqrt(sum((els-out_els_rev).^2,2));

% plot on surface to check
% figure
% plot3(els(:,1),els(:,2),els(:,3),'b.','MarkerSize',20);
% hold on;
% plot3(coord_snapped(:,1),coord_snapped(:,2),coord_snapped(:,3),'r.','MarkerSize',20);
% plot3(gs(:,1),gs(:,2),gs(:,3),'k.','MarkerSize',1);

% subfunction: Euclidean distances
function distances = dist_hermes(els, base)
distances = zeros(size(els,1),1);
for e = 1:size(els,1)
  distances(e) = sqrt(sum((els(e,:)-base').^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GridDim = determine_griddim(elec)
% assumes intact grids/strips and elec count starting at 1
% A. Stolk, 2017

% extract numbers from elec labels
digits = regexp(elec.label, '\d+', 'match');
for l=1:numel(digits)
  labels{l,1} = digits{l}{1}; % use first found digit
end

% determine grid dimensions (1st dim: number of arrays, 2nd dim: number of elecs in an array)
if isequal(numel(labels), 256)
  GridDim(1) = 16; GridDim(2) = 16;
elseif isequal(numel(labels), 64)
  GridDim(1) = 8; GridDim(2) = 8;
elseif isequal(numel(labels), 48)
  e6 = elec.elecpos(match_str(labels, num2str(6)),:);
  e7 = elec.elecpos(match_str(labels, num2str(7)),:);
  e8 = elec.elecpos(match_str(labels, num2str(8)),:);
  e9 = elec.elecpos(match_str(labels, num2str(9)),:);
  d6to7 = sqrt(sum((e6-e7).^2)); % distance of elec 6 to 7
  d8to9 = sqrt(sum((e8-e9).^2)); % distance of elec 8 to 9
  if d8to9 >= d6to7  % break between e8 and e9
    GridDim(1) = 6; GridDim(2) = 8;
  elseif d6to7 > d8to9 % break between e6 and e7
    GridDim(1) = 8; GridDim(2) = 6;
  end
elseif isequal(numel(labels), 32)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e8 = elec.elecpos(match_str(labels, num2str(8)),:);
  e9 = elec.elecpos(match_str(labels, num2str(9)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d8to9 = sqrt(sum((e8-e9).^2)); % distance of elec 8 to 9
  if d8to9 >= d4to5  % break between e8 and e9
    GridDim(1) = 4; GridDim(2) = 8;
  elseif d4to5 > d8to9 % break between e4 and e5
    GridDim(1) = 8; GridDim(2) = 4;
  end
elseif isequal(numel(labels), 24)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e6 = elec.elecpos(match_str(labels, num2str(6)),:);
  e7 = elec.elecpos(match_str(labels, num2str(7)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d6to7 = sqrt(sum((e6-e7).^2)); % distance of elec 6 to 7
  if d6to7 >= d4to5 % break between e6 and e7
    GridDim(1) = 4; GridDim(2) = 6;
  elseif d4to5 > d6to7 % break between e4 and e5
    GridDim(1) = 6; GridDim(2) = 4;
  end
elseif isequal(numel(labels), 20)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e6 = elec.elecpos(match_str(labels, num2str(6)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d5to6 = sqrt(sum((e5-e6).^2)); % distance of elec 5 to 6
  if d5to6 >= d4to5 % break between e5 and e6
    GridDim(1) = 4; GridDim(2) = 5;
  elseif d4to5 > d5to6 % break between e4 and e5
    GridDim(1) = 5; GridDim(2) = 4;
  end
elseif isequal(numel(labels), 16)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e8 = elec.elecpos(match_str(labels, num2str(8)),:);
  e9 = elec.elecpos(match_str(labels, num2str(9)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d8to9 = sqrt(sum((e8-e9).^2)); % distance of elec 8 to 9
  d4to8 = sqrt(sum((e4-e8).^2)); % distance of elec 4 to 8
  if d8to9 > 2*d4to5 % break between e8 and e9
    GridDim(1) = 2; GridDim(2) = 8;
  elseif d4to5 > 2*d4to8 % break between e4 and e5
    GridDim(1) = 4; GridDim(2) = 4;
  else
    GridDim(1) = 1; GridDim(2) = 16;
  end
elseif isequal(numel(labels), 12)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e6 = elec.elecpos(match_str(labels, num2str(6)),:);
  e7 = elec.elecpos(match_str(labels, num2str(7)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d5to6 = sqrt(sum((e5-e6).^2)); % distance of elec 5 to 6
  d6to7 = sqrt(sum((e6-e7).^2)); % distance of elec 6 to 7
  if d4to5 > 2*d5to6 % break between e4 and e5
    GridDim(1) = 3; GridDim(2) = 4; % 4x3 unsuppported
  elseif d6to7 > 2*d5to6 % break between e6 and e7
    GridDim(1) = 2; GridDim(2) = 6;
  else
    GridDim(1) = 1; GridDim(2) = 12;
  end
elseif isequal(numel(labels), 8)
  e3 = elec.elecpos(match_str(labels, num2str(3)),:);
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  d3to4 = sqrt(sum((e3-e4).^2)); % distance of elec 3 to 4
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  if d4to5 > 2*d3to4 % break between e4 and e5
    GridDim(1) = 2; GridDim(2) = 4;
  else
    GridDim(1) = 1; GridDim(2) = 8;
  end
else
  GridDim(1) = 1; GridDim(2) = numel(labels);
end

% provide feedback of what grid dimensions were found
if any(GridDim==1) % if not because of strips, this could happen in case of missing electrodes
  warning('assuming %d x %d grid dimensions: if incorrect, use cfg.pairmethod = ''pos'' instead\n', GridDim(1), GridDim(2));
else
  fprintf('assuming %d x %d grid dimensions\n', GridDim(1), GridDim(2));
end
