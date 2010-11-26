function vol = ft_headmodel_concentricspheres(geom, varargin)

% FT_HEADMODEL_CONCENTRICSPHERES
%
% Use as
%   vol = ft_headmodel_concentricspheres(geom, ...)

% get the optional arguments
conductivity = keyval('conductivity', varargin);
fitind       = keyval('fitind',       varargin); if isempty(fitind), fitind = 'all'; end

if isequal(fitind, 'all')
  fitind = 1:numel(geom);
end

% concatenate the vertices of all surfaces
pnt = [];
for i = fitind
  pnt = [pnt ; geom(i).pnt];
end

% remove double vertices
pnt  = unique(pnt, 'rows');
npnt = size(pnt, 1);

% fit a single sphere to all combined headshape points
[single_o, single_r] = fitsphere(pnt);
fprintf('initial sphere: number of unique surface points = %d\n', npnt);
fprintf('initial sphere: center = [%.1f %.1f %.1f]\n', single_o(1), single_o(2), single_o(3));
fprintf('initial sphere: radius = %.1f\n', single_r);

% fit the radius of each concentric sphere to the corresponding surface points
vol = [];
for i = 1:numel(geom)
  npnt     = size(geom(i).pnt,1);
  dist     = sqrt(sum(((geom(i).pnt - repmat(single_o, npnt, 1)).^2), 2));
  vol.r(i) = mean(dist);
end

% specify the center of the spheres
vol.o    = single_o;
vol.c    = conductivity;
vol.type = 'concentric';

% sort the spheres from the smallest to the largest
[vol.r, indx] = sort(vol.r);
vol.c = vol.c(indx);

for i=1:numel(geom)
  fprintf('concentric sphere %d: radius = %.1f, conductivity = %f\n', i, vol.r(i), vol.c(i));
end

