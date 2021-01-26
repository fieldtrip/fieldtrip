function test_issue1335

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_prepare_layout

% it is required to change to the data directory, otherwise it cannot
% automatically find the optodetemplates.xml file
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1335/'));

oxy3file = 'LR-01-2015-06-01-0002.oxy3';

matfile = 'opto.mat';
load(matfile)

ft_plot_sens(opto, 'opto', true, 'facecolor', 'k')
ft_plot_sens(opto, 'opto', false, 'facecolor', 'g')
ft_plot_axes(opto)

%%
% with the optodes from the oxy3 file, these used to be 2D

cfg = [];
cfg.opto = oxy3file;
% cfg.width = 1;
% cfg.height = 1;
% cfg.rotate = 45;

figure
ft_plot_layout(ft_prepare_layout(cfg));

%%

ft_layoutplot(cfg);

%%
% with the optodes from the opto structure, these are 3D

cfg = [];
cfg.opto = opto;
cfg.width = [];
cfg.height = [];
cfg.rotate = 90;

figure
ft_plot_layout(ft_prepare_layout(cfg));

%%

ft_layoutplot(cfg);
