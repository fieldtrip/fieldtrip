function test_issue1238

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_singleplotER ft_multiplotER ft_databrowser

% The data contains a part of a multichannel fNIRS dataset of one subject, timelocked
% for a fingertapping task (from second 0 tot 20) with in between a resting tasks.
% Before timelocking, a bandpass filter (0.01-0.14 Hz) was applied. Channels with a
% scalpcouplingindex below 0.75 were removed from the dataset (Rx4-Tx3).

% The layout for the O2Hb channels and HHb channels are plotted over each other. The
% layout was covering the left hemisphere with the outer right optodes laying over
% the midline.

% this contains layout and timelock_fingertapping
load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1238.mat'));

%% ft_singleplotER
cfg = [];
cfg.linecolor = 'rb';
cfg.colorgroups = 'sequential';

cfg.linecolor = 'r';
cfg.channel = 'Rx7-Tx7 [O2Hb]';
subplot(2,1,1); ft_singleplotER(cfg, timelock_fingertapping);

cfg.linecolor = 'b';
cfg.channel = 'Rx7-Tx7 [HHb]';
subplot(2,1,2); ft_singleplotER(cfg, timelock_fingertapping);

%% like ft_singleplotER, but multiple channels on top of each other
cfg = [];
cfg.channel = 'all';
cfg.linecolor = 'rb';
cfg.colorgroups = 'sequential';
cfg.channel = 'Rx7-Tx7*';
figure; ft_multiplotER(cfg, timelock_fingertapping);

cfg.channel = 'all';
figure; ft_multiplotER(cfg, timelock_fingertapping);

%% ft_multiplotER
cfg = [];
cfg.showlabels = 'yes';
cfg.layout = layout;
cfg.linecolor = 'rb';
cfg.colorgroups = [contains(timelock_fingertapping.label, '[O2Hb]')+2*contains(timelock_fingertapping.label, '[HHb]')];
ft_multiplotER(cfg, timelock_fingertapping);

%% ft_databrowser
cfg = [];
cfg.continuous = 'no';
cfg.linecolor = 'rb';
cfg.viewmode = 'vertical';
cfg.colorgroups = [contains(timelock_fingertapping.label, '[O2Hb]')+2*contains(timelock_fingertapping.label, '[HHb]')];
ft_databrowser(cfg, timelock_fingertapping);

