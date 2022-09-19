function failed_bug2359

% MEM 2gb
% WALLTIME 00:30:00
% DEPENDENCY ft_prepare_mesh ft_prepare_sourcemodel

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2359'));

cortex = ft_read_headshape('cortex_20484.surf.gii');
iskull = ft_read_headshape('iskull_2562.surf.gii');
oskull = ft_read_headshape('oskull_2562.surf.gii');
scalp  = ft_read_headshape('scalp_2562.surf.gii');

figure
ft_plot_mesh(cortex, 'facecolor', 'b');
ft_plot_mesh(iskull, 'facecolor', [0.5 0.5 0.5], 'edgecolor', 'none', 'facealpha', 0.5);

%% PART 1, implement the correction of the cortical sheet by moving points that are too close to the surface inward

% usually the "moveinward" cortex modification would be done for EEG-BEM, but it is
% easier to test it for a MEG singleshell model (for which it is actually not
% needed).

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, iskull);

cfg = [];
cfg.headmodel = vol;
cfg.sourcemodel = cortex;    % this is in mm
cfg.inwardshift = 0;  % this should be expressed in the units consistent with cfg.unit
cfg.moveinward = 0;
gridorig = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.headmodel = vol;
cfg.sourcemodel = cortex;      % this is in mm
cfg.inwardshift = -5;  % outward shifted
cfg.moveinward = 0;
gridoutward = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.headmodel = vol;
cfg.sourcemodel = cortex;    % this is in mm
cfg.inwardshift = 5; % inward shifted
cfg.moveinward = 0;
gridinward = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.headmodel = vol;
cfg.sourcemodel = cortex;    % this is in mm
cfg.inwardshift = 0;  % keep this at the original place
cfg.moveinward = 5;  % dipoles moved inwards
gridcorrect = ft_prepare_sourcemodel(cfg);

assert(isempty(gridorig.outside));
assert(isempty(gridoutward.outside));
assert(~isempty(gridinward.outside)); % this should have a few vertices outside the inward shifted surface
assert(isempty(gridcorrect.outside))

figure
ft_plot_mesh(iskull, 'facecolor', [0.5 0.5 0.5], 'edgecolor', 'none', 'facealpha', 0.5);
ft_plot_mesh(gridorig.pos(gridinward.outside,:), 'vertexcolor', 'b');
ft_plot_mesh(gridcorrect.pos(gridinward.outside,:), 'vertexcolor', 'r');

%% PART 2, implement the spherify of the cortical sheet

mesh = cat(1, iskull, oskull, scalp);

cfg = [];
cfg.method = 'concentricspheres';
vol = ft_prepare_headmodel(cfg, mesh);


cfg = [];
cfg.headmodel = vol;
cfg.sourcemodel = cortex;    % this is in mm
cfg.spherify = 'yes';
gridsphere = ft_prepare_sourcemodel(cfg);

assert(isempty(gridsphere.outside));

figure
ft_plot_headmodel(vol, 'edgecolor', 'none', 'facecolor', 'skin', 'facealpha', 0.5);
ft_plot_mesh(gridsphere)

%% this is a weird modification, I am just curious to see how it works

cfg = [];
cfg.method = 'singlesphere';
vol = ft_prepare_headmodel(cfg, iskull);

cfg = [];
cfg.headmodel = vol;
cfg.sourcemodel.xgrid = -200:10:200;    % this is in mm
cfg.sourcemodel.ygrid = -200:10:200;    % this is in mm
cfg.sourcemodel.zgrid =  -50:10:150;    % this is in mm
cfg.spherify = 'yes';
gridsphere = ft_prepare_sourcemodel(cfg);

figure
ft_plot_headmodel(vol, 'edgecolor', 'none', 'facecolor', 'skin', 'facealpha', 0.5);
ft_plot_mesh(gridsphere)


