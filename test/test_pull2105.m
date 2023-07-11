function test_pull2105

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY ft_topoplotER ft_plot_topo
% NODATA

% create artificial data with 2 channels
data.time = [-10: 10];
data.label = {'C3'; 'C4'};
data.avg = [ones(1,21)*1; ones(1,21)*-1];
data.dimord = 'chan_time';

% create layout with two masks
layout.pos = [-1 0; +1 0];
layout.label = {'C3'; 'C4'};
layout.width = [0.2; 0.2];
layout.height = [0.2; 0.2];
layout.mask = {[-1.5 -0.5; -0.5 -0.5; -0.5 0.5; -1.5 0.5; -1.5 -0.5] [0.5 -0.5; 1.5 -0.5; 1.5 0.5; 0.5 0.5; 0.5 -0.5;]};
layout.outline = {[-2 -1; 2 -1; 2 1; -2 1]};

% show layout
figure; ft_plot_layout(layout)

% first topoplot
cfg =[];
cfg.xlim = [1 5];
cfg.zlim = [-1 1];
cfg.layout = layout;
cfg.interplimits = 'mask_individual';
ft_topoplotER(cfg, data)
% This creates a figure that contains a blue and a yellow mask (as it
% should)


% second topoplot
cfg =[];
cfg.xlim = [1 5];
cfg.zlim = [-1 1];
cfg.layout = layout;
cfg.interplimits = 'mask_individual';
ft_topoplotER(cfg, data)
% This create a figure with an interpolation of the data between
% the mask (which should not be the case).

