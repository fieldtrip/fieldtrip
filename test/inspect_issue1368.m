function inspect_issue1368

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_plot_sens

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1368/hull.mat'))
load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1368/elec.mat'))

ft_plot_sens(elec, 'elecshape', 'disc', 'headshape', hull);