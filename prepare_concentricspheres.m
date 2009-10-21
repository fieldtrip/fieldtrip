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
% $Log: prepare_concentricspheres.m,v $
% Revision 1.8  2009/07/16 09:14:52  crimic
% part of code reimplemented in function prepare_mesh_headshape.m
%
% Revision 1.7  2009/06/23 14:59:28  crimic
% use of plotting toolbox funtion: plot_mesh
%
% Revision 1.6  2009/05/29 11:40:07  roboos
% only convert cfg.headshape from config to struct in case it is present
%
% Revision 1.5  2009/05/25 08:05:18  roboos
% ensure that cfg.headshape is a sturct and not a config object (in case tracking is on)
%
% Revision 1.4  2009/05/14 19:21:36  roboos
% consistent handling of cfg.headshape in code and documentation
%
% Revision 1.3  2009/04/01 12:28:58  roboos
% use Taubin's method instead of nonlinear search (thanks to Jean and Guillaume)
%
% Revision 1.2  2009/02/05 10:22:55  roboos
% don't open new figure, clear the existing one
%
% Revision 1.1  2009/01/05 13:06:39  roboos
% initial version of Vladimir with some extensions/improvements
%

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

