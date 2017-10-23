function test_bug1651

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_prepare_mesh ft_datatype_segmentation

readfromdisk = true;

if readfromdisk
  % the segmentation takes quite some time
  % furthermore we want to ensure that we have the old-stype segmentation
  load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1651.mat'));
else
  % this is the original code to create the segmentation
  % create a 3 layered segmentation
  mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));
  
  % this speeds up the subsequent stuff
  cfg = [];
  cfg.downsample = 2;
  mri = ft_volumedownsample(cfg, mri);
  
  cfg = [];
  cfg.output = {'brain', 'scalp', 'skull'};
  seg2 = ft_volumesegment(cfg, mri);
end

% the following results in the error
cfg = [];
cfg.tissue = {'brain', 'skull', 'scalp'};
cfg.numvertices = [3000 2000 1000];
bnd = ft_prepare_mesh(cfg, seg2);

figure
ft_plot_mesh(bnd)
alpha 0.1
