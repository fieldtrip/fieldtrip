function vol = ft_headmodel_concentricspheres(geom, varargin)

% FT_HEADMODEL_CONCENTRICSPHERES creates a volume conduction model
% of the head based on three or four concentric spheres. For a 3-sphere
% model the spheres represent the skin surface, the outside of the
% skull and the inside of the skull For a 4-sphere model, the surfaces
% describe the skin, the outside-skull, the inside-skull and the inside of the 
% cerebro-spinal fluid (CSF) boundaries.
% 
% The innermost surface is sometimes also referred to as the brain
% surface, i.e. as the outside of the brain volume.
% 
% This function takes as input a single headshape described with
% points and fits the spheres to this surface. If you have a set of
% points describing each surface, then this function fits the spheres
% to all individual surfaces.
% 
% Use as
%   vol = ft_headmodel_concentricspheres(geom, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   conductivity     = vector with the conductivity of each compartment
%   fitind           = vector with indices of the surfaces to use in fitting the center of the spheres
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

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

