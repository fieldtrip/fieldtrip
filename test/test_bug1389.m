function test_bug1389

% MEM 1500mb
% WALLTIME 00:10:00

% there is confusion as to at what level in the cfg tre preproc options
% need to be specified: for ft_preprocessing this is to be at the main
% level of the cfg, for all other functions that rely on preproc, it should
% be in cfg.preproc.XXX

ok = true;

data = [];
data.trial{1} = randn(2,100)+10;
data.time{1}  = (-20:79)./100;
data.label    = {'1';'2'};

% this is how it should work
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-inf 0];
data2 = ft_preprocessing(cfg, data);
ok = ok && all(mean(data.trial{1},2)~=mean(data2.trial{1},2));

% this is how it shouldn't work
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindwo = [-inf 0];
data2b = ft_preprocessing(cfg, data);
ok = ok && all(mean(data.trial{1},2)==mean(data2b.trial{1},2));

% this is how it should work
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-inf 0];
data3 = ft_timelockanalysis(cfg, data);
ok = ok && all(mean(data.trial{1},2)~=mean(data3.avg,2));
