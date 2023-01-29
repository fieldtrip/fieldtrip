function test_bug2990

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_prepare_sourcemodel ft_volumereslice

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2990/4roboos_part1.mat'));
% the mri is from the dicom files, mri_realigned is in CTF coordinates, mri_realigned_resliced is after reslicing

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2990/4roboos_part2.mat'));
% the cfg contains cfg.mri, which is more or less the mri_realigned (but more accurately coregistered than the sloppy ones in part1)
% the cfg also contains cfg.grid

cfg.checkconfig = 'loose';

cfg1 = cfg;
cfg1.mri = mri_realigned;
grid1 = ft_prepare_sourcemodel(cfg1);

cfg2 = cfg;
cfg2.mri = mri_realigned_resliced;
grid2 = ft_prepare_sourcemodel(cfg2);

figure; ft_plot_mesh(cfg.grid.template.pos(cfg.grid.template.inside,:))
figure; ft_plot_mesh(grid1.pos(grid1.inside,:))
figure; ft_plot_mesh(grid2.pos(grid2.inside,:))

% cfg.sourcemodel.template.pos(1,:)
% ans =
%     -8   -11    -7
%
% grid1.pos(1,:)
% ans =
%    -5.3781    8.3252   -0.2330
%
% Hmm, the first is in MNI, the second in CTF coordinates. Both are left-anterior-inferior.
