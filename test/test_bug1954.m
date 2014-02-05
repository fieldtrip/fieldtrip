function test_bug1954

% MEM 2gb
% WALLTIME 00:10:00

load /home/common/matlab/fieldtrip/template/headmodel/standard_mri.mat


cfg           = [];
cfg.output    = {'brain','skull','scalp'};
cfg.coordsys  = 'ras';
segmentedmri  = ft_volumesegment(cfg, mri);
% I get warning (exact file name changes each time):
% Warning: could not open /tmp/tpe20fae07_8585_42ab_a6e4_0e4c656f97b3.img 

if 0 % so this doesn't run during automated testing
  figure;imagesc(squeeze(segmentedmri.scalp(:,110,:)))
  % voxels on top of head included in scalp
end

cfg=[];
cfg.numvertices=3000;
bnd=ft_prepare_mesh(cfg,segmentedmri);
if 0 
 figure;ft_plot_mesh(bnd(1),'facealpha',.2);hold on;ft_plot_mesh(bnd(2))
 figure;ft_plot_mesh(bnd(1),'facealpha',.2);hold on;ft_plot_mesh(bnd(3))
 figure;ft_plot_mesh(bnd(3),'facealpha',.2);hold on;ft_plot_mesh(bnd(2))
 % order is now 1) scalp, 2) brain, 3) skull
end


cfg = [];
cfg.method = 'bemcp';
vol1 = ft_prepare_headmodel(cfg, segmentedmri);
% vol1.mat are all NaN

cfg = [];
cfg.method = 'dipoli';
vol2 = ft_prepare_headmodel(cfg, segmentedmri);
% looks ok, but no .mat to assess

segmentedmri.bnd=bnd;
cfg=[];
cfg.method='bemcp';
vol1o=ft_prepare_bemmodel(cfg,segmentedmri);
% vol1o.mat also all NaN


