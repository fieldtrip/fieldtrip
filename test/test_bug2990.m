function test_bug2990

% WALLTIME 0:10:00
% MEM 2gb

% TEST ft_prepare_sourcemodel ft_volumereslice

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2990/4roboos_part1'));
% the mri is from the dicom files, mri_realigned is in CTF coordinates, mri_realigned_resliced is after reslicing

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2990/4roboos_part2'));
% the cfg contains cfg.mri, which is more or less the mri_realigned (but more accurately coregistered than the sloppy ones in part1)

cfg1 = cfg;
cfg1.mri = mri_realigned;
grid1 = ft_prepare_sourcemodel(cfg);

cfg2 = cfg;
cfg2.mri = mri_realigned_resliced;
grid2 = ft_prepare_sourcemodel(cfg);

assert(all(abs(grid1.pos(:)-grid2.pos(:))<0.001));

figure; ft_plot_mesh(cfg.grid.template.pos(cfg.grid.template.inside,:))
figure; ft_plot_mesh(grid1.pos(grid1.inside,:))
figure; ft_plot_mesh(grid2.pos(grid2.inside,:))

% cfg.grid.template.pos(1,:)
% ans =
%     -8   -11    -7
%
% grid1.pos(1,:)
% ans =
%    -5.3781    8.3252   -0.2330
%
% Hmm, the first is in MNI, the second in CTF coordinates. Both are left-anterior-inferior.
