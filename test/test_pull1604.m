function test_pull1604

% MEM 8gb
% WALLTIME 00:40:00
% DEPENDENCY ft_inside_headmodel

datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/');
load(fullfile(datadir, 'segmentedmri.mat'));
load(fullfile(datadir, 'sourcemodel.mat'));

cfg = [];
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg,segmentedmri);

headmodel      = mesh;
headmodel.type = 'simbio';
headmodel.conductivity = zeros(size(headmodel.tissue)); % this is irrelevant for now.

headmodel = ft_convert_units(headmodel, 'cm');

% compare three different ways of determining the indices of the closest
% points, my_dsearchn is a copy of the subfunction in ft_inside_headmodel:
tic; indx1 = dsearchn(headmodel.pos, sourcemodel.pos(sourcemodel.inside,:));   t(1) = toc; 
tic; indx2 = my_dsearchn(headmodel.pos, sourcemodel.pos(sourcemodel.inside,:), 0); t(2) = toc;
tic; indx2b = my_dsearchn(headmodel.pos, sourcemodel.pos(sourcemodel.inside,:)); t(3) = toc;
tic; indx3 = knnsearch(headmodel.pos, sourcemodel.pos(sourcemodel.inside,:));   t(4) = toc;
assert(isequal(indx1, indx2));
assert(isequal(indx2, indx3));
assert(isequal(indx2, indx2b));

fprintf('dsearchn took %3.2f seconds\n',t(1));
fprintf('native use of my_dsearchn took %3.2f seconds\n',t(2));
fprintf('my_dsearchn took %3.2f seconds\n', t(3));
fprintf('knnsearchn took %3.2f seconds\n',t(4));

% run ft_inside_headmodel once just to check that it runs
tic; inside = ft_inside_headmodel(sourcemodel.pos, headmodel); t(5) = toc;
fprintf('ft_inside_headmodel took %3.2f seconds\n',t(5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local subfunction that is much faster than dsearchn

function indx = my_dsearchn(pos1, pos2, flag)

% indx = zeros(size(pos2,1),1);
% mind = inf(size(pos2,1),1);
% for k = 1:size(pos1,1)
%   dpos = bsxfun(@minus, pos2, pos1(k,:));
%   thisd = sum(dpos.^2,2);
%   issmaller = thisd<mind;
%   mind(issmaller) = thisd(issmaller);
%   indx(issmaller) = k;
% end

% the idea is that the distance between 2 points is:
% 
% sqrt(sum((p1(x,y,z)-p2(x,y,z)).^2)
% 
% since we are dealing with relative distances, we can get rid of the sqrt:
% so we need to compute: 
%
% sum((p1(x,y,z)-p2(x,y,z)).^2)
%
% this is the same as:
%
% (p1x-p2x)^2 + (p1y-p2y)^2 + (p1z-p2z)^2
%
% or, equivalently:
%
% p1x^2 + p2x^2 - 2*p1x*p2x+ ... 
% 
% reordering:
%
% (p1x^2 + p1y^2 + p1z^2) + cross-terms + (p2x^2 + p2y^2 + p2z^2)
%
% the last term between brackets is the same for each position-of-interest:
% so it does not change the relative distance, and the first term between
% brackets only needs to be computed once (below denoted as the 'offset'
% variable.

if nargin<3
  flag = true;
end

if flag && exist('knnsearch', 'file')
  % use much faster knnsearch if available on the path
  indx = knnsearch(pos1, pos2);
  return;
end

offset = (pos1.^2)*[1;1;1];

% not sure whether this speeds up things, but it does not hurt do the
% operations in order of overall offset (i.e. absolute distance of the
% headmodel points to the origin)
[srt, ix] = sort(-offset);
indx   = zeros(size(pos2,1),1);
mind   = inf(1,size(pos2,1));

% transpose once, to speed up matrix computations
pos2 = pos2';

chunksize = 250; % just a number, could be optimized
chunks    = [(0:chunksize:(numel(offset)-1)) numel(offset)];

% loop across blocks of headmodel points, and iteratively update the
% index to the nearest sourcemodel point, based on the shortcut heuristic
% explained above

for k = 1:(numel(chunks)-1)
  iy = ix((chunks(k)+1):chunks(k+1));
  
  %thisd = offset(k) - 2.*(pos2(:,1).*pos1(k,1)+pos2(:,2).*pos1(k,2)+pos2(:,3).*pos1(k,3));
  thisd  = offset(iy)./2 - pos1(iy,:)*pos2;%pos2(:,1)*pos1(k,1)-pos2(:,2)*pos1(k,2)-pos2(:,3)*pos1(k,3);
  [m, i] = min(thisd, [], 1);
  issmaller = m<mind;
  mind(issmaller) = m(issmaller);
  indx(issmaller) = iy(i(issmaller));
end

