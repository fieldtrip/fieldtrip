function test_issue714

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_databrowser ft_multiplotER ft_singleplotER

% hand-craft a dataset with two trials/conditions and 6 channels
% the six channels come in three a/b pairs, resembing NIRS oxi/deoxihemoglobin channels
data = [];
data.label = {
  '1a'
  '2a'
  '3a'
  '1b'
  '2b'
  '3b'};
data.chantype = {
  'oxi'
  'oxi'
  'oxi'
  'deoxi'
  'deoxi'
  'deoxi'};

data.trialinfo = [1 2]'; % for selective averaging
data.time{1} = (1:400)/400-0.1;
data.time{2} = (1:400)/400-0.1;

data.trial{1} = zeros(6,400);
data.trial{1}(1, 51:150) =  1;
data.trial{1}(2,151:250) =  1;
data.trial{1}(3,251:350) =  1;
data.trial{1}(4, 51:150) = -1;
data.trial{1}(5,151:250) = -1;
data.trial{1}(6,251:350) = -1;
data.trial{1} = data.trial{1} + randn(6,400)*0.02;

data.trial{2} = zeros(6,400);
data.trial{2}(1, 51:150) =  1.2;
data.trial{2}(2,151:250) =  1.2;
data.trial{2}(3,251:350) =  1.2;
data.trial{2}(4, 51:150) = -0.8;
data.trial{2}(5,151:250) = -0.8;
data.trial{2}(6,251:350) = -0.8;
data.trial{2} = data.trial{2} + randn(6,400)*0.02;

cfg = [];
cfg.trials = data.trialinfo==1;
cond1 = ft_timelockanalysis(cfg, data);
cfg.trials = data.trialinfo==2;
cond2 = ft_timelockanalysis(cfg, data);

cond1.chantype = data.chantype;
cond2.chantype = data.chantype;


%%

% make a simple layout with the "a" channels overlapping the "b" channels
% the outline and mask are not updated, so it is not so nice yet

cfg = [];
cfg.layout = 'vertical';
cfg.skipscale = 'no';
cfg.skipcomnt = 'no';
layout = ft_prepare_layout(cfg, data);
layout.pos(4:6,:) = layout.pos(1:3,:);

ft_plot_layout(layout)

%%
% this uses the old options

cfg = [];
cfg.checkconfig = 'loose'; % this uses an old option
cfg.viewmode = 'vertical';
cfg.channelcolormap = colormap('lines');

cfg.colorgroups = 'sequential';
ft_databrowser(cfg, data);

cfg.colorgroups = 'allblack';
ft_databrowser(cfg, data);

cfg.colorgroups = 'labelchar2';
ft_databrowser(cfg, data);

cfg.colorgroups = 'chantype';
ft_databrowser(cfg, data);

cfg.colorgroups = [1 1 1 2 2 2 ];
ft_databrowser(cfg, data);

%%
% this uses the new options

cfg = [];
cfg.layout = layout;
% cfg.viewmode = 'butterfly';
cfg.viewmode = 'vertical';
cfg.colorgroups = 'sequential';

cfg.linewidth = 2;
cfg.linestyle = '-';

% cfg.linecolor = colormap('lines');
ft_databrowser(cfg, data);

cfg.linecolor = 'brgymc';
ft_databrowser(cfg, data);

cfg.linecolor = 'm';
ft_databrowser(cfg, data);

%%

cfg = [];
cfg.checkconfig = 'loose'; % this uses an old option
cfg.linestyle = '-';
cfg.linewidth = 0.5;
cfg.graphcolor = 'mc';
cfg.channel = 1;
figure
ft_singleplotER(cfg, cond1, cond2);

%%
% this uses the new options AND some new functionality

cfg = [];
cfg.linestyle = {'-'};
cfg.linewidth = 2;
cfg.linecolor = 'mc';
cfg.showlegend = 'yes';
cfg.channel = 1;
figure
ft_singleplotER(cfg, cond1, cond2);

cfg.maskstyle = 'difference';
figure
ft_singleplotER(cfg, cond1, cond2);


%%

cfg = [];
cfg.checkconfig = 'loose'; % this uses an old option
cfg.layout = layout;
cfg.linestyle = '-';
cfg.linewidth = 0.5;
cfg.graphcolor = 'mc';
ft_multiplotER(cfg, cond1, cond2);

%%
% this uses the new options

cfg = [];
cfg.layout = layout;
cfg.linestyle = {'-'};
cfg.linewidth = 2;
cfg.linecolor = 'mc';
cfg.showlegend = 'yes';
figure
ft_multiplotER(cfg, cond1, cond2);

cfg.maskstyle = 'difference';
figure
ft_multiplotER(cfg, cond1, cond2);

%%
% this uses the new options AND some new functionality

cfg = [];
cfg.layout = layout;
cfg.linecolor = 'rgbymc';
cfg.colorgroups = 'condition';
figure; ft_multiplotER(cfg, cond1);
cfg.colorgroups = 'sequential';
figure; ft_multiplotER(cfg, cond1);
cfg.colorgroups = 'allblack';
figure; ft_multiplotER(cfg, cond1);
cfg.colorgroups = 'labelchar2';
figure; ft_multiplotER(cfg, cond1);
cfg.colorgroups = 'chantype';
figure; ft_multiplotER(cfg, cond1);

%%
% this combines new (color) with old (bufferfly) functionality

cfg = [];
cfg.layout = 'butterfly';
cfg.linecolor = 'rgbymc';
cfg.colorgroups = 'sequential';
figure; ft_multiplotER(cfg, cond1);

