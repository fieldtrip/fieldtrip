function test_bug1409

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_selectdata

% the issue is that ft_selectdata_new removes the dimord

data = [];
data.trial{1} = randn(2,100);
data.time{1}  = [0:99]./100;
data.label    = {'chan1';'chan2'};

data = ft_checkdata(data, 'datatype', 'timelock');

cfg = [];
cfg.channel = 'chan1';
data2 = ft_selectdata(cfg, data);

assert(isfield(data2,'dimord'));
