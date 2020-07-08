function test_ft_appendlayout

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

nchan = 24;
ntime = 1000;
dataA = [];
for i=1:nchan
  dataA.label{i} = sprintf('A%02d', i);
end
dataA.avg = rand(nchan, ntime);
dataA.time = (1:ntime)/ntime;

nchan = 16;
ntime = 1000;
dataB = [];
for i=1:nchan
  dataB.label{i} = sprintf('B%02d', i);
end
dataB.avg = rand(nchan, ntime);
dataB.time = (1:ntime)/ntime;

nchan = 8;
ntime = 1000;
dataC = [];
for i=1:nchan
  dataC.label{i} = sprintf('C%02d', i);
end
dataC.avg = rand(nchan, ntime);
dataC.time = (1:ntime)/ntime;

nchan = 10;
ntime = 1000;
dataD = [];
for i=1:nchan
  dataD.label{i} = sprintf('D%02d', i);
end
dataD.avg = rand(nchan, ntime);
dataD.time = (1:ntime)/ntime;

% append the iEEG data from two grids and two shafts
cfg = [];
data = ft_appenddata(cfg, dataA, dataB, dataC, dataD);

% compute the average
cfg = [];
timelock = ft_timelockanalysis(cfg, data);

%%

cfg = [];
cfg.channel = 'A*';
cfg.layout = 'ordered';
cfg.rows = 3;
cfg.skipscale = 'no';
cfg.skipcomnt = 'no';
layoutA = ft_prepare_layout(cfg, dataA);
ft_layoutplot(cfg, dataA);

%%

cfg = [];
cfg.channel = 'B*';
cfg.layout = 'ordered';
cfg.direction = 'BTRL';
cfg.rows = 4;
cfg.skipscale = 'no';
cfg.skipcomnt = 'no';
layoutB = ft_prepare_layout(cfg, dataB);
ft_layoutplot(cfg, dataB);

%%

cfg = [];
cfg.channel = 'C*';
cfg.layout = 'vertical';
cfg.skipscale = 'no';
cfg.skipcomnt = 'no';
layoutC = ft_prepare_layout(cfg, dataC);
ft_layoutplot(cfg, dataC);

%%

cfg = [];
cfg.channel = 'D*';
cfg.layout = 'vertical';
cfg.direction = 'BT';
cfg.skipscale = 'no';
cfg.skipcomnt = 'no';
layoutD = ft_prepare_layout(cfg, dataD);
ft_layoutplot(cfg, dataD);

%%

cfg = []; cfg.direction = 'horizontal'; cfg.align = 'center';
combined = ft_appendlayout(cfg, layoutA, layoutB, layoutC, layoutD);
figure; ft_plot_layout(combined); title([cfg.direction ' ' cfg.align]);

cfg = [];
cfg.layout = combined;
cfg.commentpos = 'lefttop';
cfg.scalepos = 'righttop';
ft_multiplotER(cfg, timelock);

%%

cfg = []; cfg.direction = 'horizontal'; cfg.align = 'top';
combined = ft_appendlayout(cfg, layoutA, layoutB, layoutC, layoutD);
figure; ft_plot_layout(combined); title([cfg.direction ' ' cfg.align]);

%%

cfg = []; cfg.direction = 'horizontal'; cfg.align = 'bottom';
combined = ft_appendlayout(cfg, layoutA, layoutB, layoutC, layoutD);
figure; ft_plot_layout(combined); title([cfg.direction ' ' cfg.align]);

%%

cfg = []; cfg.direction = 'vertical'; cfg.align = 'center';
combined = ft_appendlayout(cfg, layoutA, layoutB, layoutC, layoutD);
figure; ft_plot_layout(combined); title([cfg.direction ' ' cfg.align]);

%%

cfg = []; cfg.direction = 'vertical'; cfg.align = 'left';
combined = ft_appendlayout(cfg, layoutA, layoutB, layoutC, layoutD);
figure; ft_plot_layout(combined); title([cfg.direction ' ' cfg.align]);

%%

cfg = []; cfg.direction = 'vertical'; cfg.align = 'right';
combined = ft_appendlayout(cfg, layoutA, layoutB, layoutC, layoutD);
figure; ft_plot_layout(combined); title([cfg.direction ' ' cfg.align]);
