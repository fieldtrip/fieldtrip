function test_bug996

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_multiplotER ft_singleplotER

% function to reproduce bug 996 and test the fix
% the problem is the following: JM removed the occurrence of
% cfg.xparam/cfg.yparam (since they are not specifiable), and changed it
% into xparam and yparam. These are defined as strings, but later on xparam
% is replaced by a vector (defining the axis of the xparam). When multiple
% inputs are present, the next iteration in the for-loop expects xparam
% again to be a string, which is not anymore the case

% create some data
cfg = [];
cfg.layout = 'CTF151.lay';
lay = ft_prepare_layout(cfg);

x.time  = 1:10;
x.label = ft_channelselection('MEG', lay.label);
x.avg   = rand(151,10);
x.dimord = 'chan_time';

y.time  = 1:10;
y.label = ft_channelselection('MEG', lay.label);
y.avg   = rand(151,10);
y.dimord = 'chan_time';

figure;ft_multiplotER(cfg,x,y);

cfg.channel = {'MLC11'};
figure;ft_singleplotER(cfg,x,y);
