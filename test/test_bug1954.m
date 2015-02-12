function test_bug1954

% MEM 2gb
% WALLTIME 00:10:00

load(fullfile(fileparts(which('ft_defaults')), 'template', 'headmodel', 'standard_mri.mat'));

mri.coordsys = 'ras'; % this can also be determined with ft_determine_coordsys

cfg = [];
cfg.output = {'brain','skull','scalp'};
% cfg.coordsys = 'ras'; % not supported any more, should be specified in the input data
segmentedmri = ft_volumesegment(cfg, mri);

% I get warning (exact file name changes each time):
% Warning: could not open /tmp/tpe20fae07_8585_42ab_a6e4_0e4c656f97b3.img

figure; imagesc(squeeze(segmentedmri.scalp(:,110,:)))
% voxels on top of head included in scalp

cfg = [];
cfg.numvertices = 1000;
cfg.tissue = {'scalp', 'skull', 'brain'};
bnd = ft_prepare_mesh(cfg,segmentedmri);

figure;ft_plot_mesh(bnd(1),'facecolor', 'none'); hold on; ft_plot_mesh(bnd(2),'facecolor', 'none')
figure;ft_plot_mesh(bnd(1),'facecolor', 'none'); hold on; ft_plot_mesh(bnd(3),'facecolor', 'none')
figure;ft_plot_mesh(bnd(2),'facecolor', 'none'); hold on; ft_plot_mesh(bnd(3),'facecolor', 'none')
% order is now 1) scalp, 2) brain, 3) skull

cfg = [];
cfg.method = 'bemcp';
cfg.numvertices = 1000;
cfg.tissue = {'scalp', 'skull', 'brain'};
vol1 = ft_prepare_headmodel(cfg, segmentedmri);
assert(~any(isnan(vol1.mat(:))), 'there is a NaN in vol1.mat');

cfg = [];
cfg.method = 'dipoli';
vol2 = ft_prepare_headmodel(cfg, segmentedmri);
assert(~any(isnan(vol2.mat(:))), 'there is a NaN in vol2.mat');

cfg = [];
cfg.method = 'bemcp';
vol3 = ft_prepare_headmodel(cfg,bnd);
any(isnan(vol3.mat(:)))
assert(~any(isnan(vol3.mat(:))), 'there is a NaN in vol3.mat');
