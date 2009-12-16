function [vol, cfg] = prepare_concentricspheres(cfg)

% PREPARE_CONCENTRICSPHERES creates a EEG volume conductor model with
% multiple concentric spheres.
%
% Use as
%   [vol, cfg] = prepare_concentricspheres(cfg)
%
% The input configuration should contain
%   cfg.headshape     = a filename containing headshape, a Nx3 matrix with surface 
%                       points, or a structure with a single or multiple boundaries
%   cfg.conductivity  = conductivity values for the model (default = [0.3300 1 0.0042 0.3300])
%   cfg.fitind        = indices of shapes to use for fitting the center (default = 'all')
%   cfg.nonlinear     = 'yes' or 'no' (default = 'yes')
%   cfg.feedback      = 'yes' or 'no' (default = 'yes')
%
% Example:
%
%   % first create 4 surfaces that represent the brain, csf, skull and skin
%   radius = [86 88 92 100];
%   headshape = [];
%   for i=1:4
%     pnt = randn(100,3);
%     for j=1:size(pnt,1)
%       pnt(j,:) = pnt(j,:) ./ norm(pnt(j,:));
%     end
%     headshape(i).pnt = radius(i) .* pnt + 0.1*randn(size(pnt));
%   end
%
%   % then construct a volume conduction model of the head by fitting 4 concentric spheres
%   cfg = [];
%   cfg.headshape    = headshape;
%   cfg.conductivity = [0.3300 1 0.0042 0.3300]
%   [vol, cfg] = prepare_concentricspheres(cfg)

% Copyright (C) 2009, Vladimir Litvak & Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'fitind'),        cfg.fitind = 'all';                            end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'yes';                          end
if ~isfield(cfg, 'conductivity'),  cfg.conductivity = [0.3300 1 0.0042 0.3300];   end
if ~isfield(cfg, 'numvertices'),   cfg.numvertices = 'same';                      end

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested config-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

cfg = checkconfig(cfg, 'forbidden', 'nonlinear');

% get the surface describing the head shape
headshape = prepare_mesh_headshape(cfg);

if strcmp(cfg.fitind, 'all')
  fitind = 1:numel(headshape);
else
  fitind = cfg.fitind;
end

% concatenate the vertices of all surfaces
pnt = [];
for i = fitind
  pnt = [pnt ; headshape(i).pnt];
end
% remove double vertices
pnt = unique(pnt, 'rows');

Npnt = size(pnt, 1);

% set up an empty figure
if strcmp(cfg.feedback, 'yes')
  clf
  hold on
  axis equal
  axis vis3d
  axis off
  drawnow
  colors = {'b', 'y', 'm', 'r'};
  [sphere_pnt, sphere_tri] = icosahedron162;
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(pnt);
fprintf('initial sphere: number of surface points = %d\n', Npnt);
fprintf('initial sphere: center = [%.1f %.1f %.1f]\n', single_o(1), single_o(2), single_o(3));
fprintf('initial sphere: radius = %.1f\n', single_r);

% fit the radius of each concentric sphere to the corresponding surface points
vol = [];
vol.o = single_o;
for i = 1:numel(headshape)
  dist     = sqrt(sum(((headshape(end-i+1).pnt - repmat(single_o, size(headshape(end-i+1).pnt,1), 1)).^2), 2));
  vol.r(i) = mean(dist);

  if strcmp(cfg.feedback, 'yes')
    if ~isfield(headshape(end-i+1), 'tri')
      headshape(end-i+1).tri = [];
    end

    % plot the original surface
    bndtmp = [];
    bndtmp.pnt = headshape(end-i+1).pnt;
    bndtmp.tri = headshape(end-i+1).tri;
    plot_mesh(bndtmp,'facecolor','none')

    % plot the sphere surface
    bndtmp = [];
    bndtmp.pnt = sphere_pnt*vol.r(i) + repmat(single_o, size(sphere_pnt, 1), 1);
    bndtmp.tri = sphere_tri;
    plot_mesh(bndtmp,'edgecolor',colors{mod(i, numel(colors)) + 1},'facecolor','none');
  end
end

if numel(cfg.conductivity)==numel(headshape)
  vol.c = cfg.conductivity;
else
  error('incorrect specification of cfg.conductivity');
end

vol.type = 'concentric';

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

