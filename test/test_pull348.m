function test_pull348

% MEM 1000mb
% WALLTIME 00:10:00
%
% TEST ft_componentanalysis 
% TEST bsscca

% this function tests whether the bsscca method works ok in
% ft_componentanalysis

% generate some data
data = [];
data.trial{1} = randn(3,100);
data.trial{2} = randn(3,100);
data.time{1}  = (0:99)./100;
data.time{2}  = (0:99)./100;
data.label    = {'chan1';'chan2';'chan3'};

refdata = [];
refdata.trial{1} = randn(1,100);
refdata.trial{2} = randn(1,100);
refdata.time     = data.time;
refdata.label    = {'refchan'};

% method 1
cfg = [];
cfg.method = 'bsscca';
cfg.cellmode = 'yes';
cfg.updatesens = 'no';
cfg.bsscca.delay = 1;
comp = ft_componentanalysis(cfg, data);
assert(ft_datatype(comp, 'comp'));

% method 2
cfg.bsscca.refdata = refdata.trial;
cfg.bsscca.refdelay = (0:10);
comp2 = ft_componentanalysis(cfg, data);
assert(ft_datatype(comp2, 'comp'));

