function test_bug1954

% MEM 2gb
% WALLTIME 00:45:00

[ftver, ftpath] = ft_version;
load(fullfile(ftpath, 'template', 'headmodel', 'standard_mri.mat'));

mri.coordsys = 'ras'; % this can also be determined with ft_determine_coordsys

cfg = [];
cfg.output = {'brain','skull','scalp'};
% cfg.coordsys = 'ras'; % not supported any more, should be specified in the input data
segmentedmri = ft_volumesegment(cfg, mri);

scalp = imfill(segmentedmri.scalp, [128 128 128]);
skull = imfill(segmentedmri.skull, [128 128 128]);
brain = imfill(segmentedmri.brain, [128 128 128]);

% the brain and skull go too far down at the bottom
skull(:,:,1:5) = 0;
brain(:,:,1:10) = 0;

% ensure that they don't overlap
skull = skull & imerode(scalp, strel_bol(2));
brain = brain & imerode(skull, strel_bol(2));

segmentedmri.scalp = scalp;
segmentedmri.skull = skull;
segmentedmri.brain = brain;

if false
  segmentedmri.seg = segmentedmri.scalp + segmentedmri.skull + segmentedmri.brain;
  cfg = [];
  cfg.funparameter = 'seg';
  ft_sourceplot(cfg, segmentedmri);
end

% I get warning (exact file name changes each time):
% Warning: could not open /tmp/tpe20fae07_8585_42ab_a6e4_0e4c656f97b3.img

figure; imagesc(squeeze(segmentedmri.scalp(:,110,:)))
% voxels on top of head included in scalp

cfg = [];
cfg.numvertices = 1000;
cfg.tissue = {'scalp', 'skull', 'brain'};
bnd = ft_prepare_mesh(cfg,segmentedmri);

figure; ft_plot_mesh(bnd(1),'facecolor', 'none'); hold on; ft_plot_mesh(bnd(2),'facecolor', 'none')
figure; ft_plot_mesh(bnd(1),'facecolor', 'none'); hold on; ft_plot_mesh(bnd(3),'facecolor', 'none')
figure; ft_plot_mesh(bnd(2),'facecolor', 'none'); hold on; ft_plot_mesh(bnd(3),'facecolor', 'none')
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
assert(~any(isnan(vol3.mat(:))), 'there is a NaN in vol3.mat');
