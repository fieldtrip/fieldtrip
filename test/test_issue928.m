function test_issue928

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_singleplotER ft_plot_vector

% create some data
data1 = [];
data1.avg  = randn(2,100);
data1.time = (0:99)./100;
data1.dimord = 'chan_time';
data1.label = {'chan01';'chan02'};

data2 = [];
data2.avg  = randn(2,100);
data2.time = (0:99)./100;
data2.dimord = 'chan_time';
data2.label = {'chan01';'chan02'};

cfg = [];
cfg.parameter = 'avg';
cfg.showlabels= 'yes';
cfg.linestyle = {'-' '--'};

ft_singleplotER(cfg, data1, data2);
h = get(gca, 'children');

assert(~isequal(get(h(1),'linestyle'), get(h(2),'linestyle')));
