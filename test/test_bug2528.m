function test_bug2528

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_math ft_selectdata ft_selectdata_new

% data1 =
%           avg: [63x2000 double]
%           var: [63x2000 double]
%          time: [1x2000 double]
%           dof: [63x2000 double]
%         label: {63x1 cell}
%         trial: [42x63x2000 double]
%        dimord: 'rpt_chan_time'
%     trialinfo: [42x5 double]
%           cfg: [1x1 struct]

data1 = [];
data1.avg = randn(63,2000);
data1.var = randn(63,2000);
data1.dof = randn(63,2000);
data1.time = (1:2000)/1000;
for i=1:63
  data1.label{i} = num2str(i);
end
data1.trial = randn(42,63,2000);
data1.trialinfo = randn(42,5);
data1.dimord = 'rpt_chan_time';

data2 = [];
data2.avg = randn(63,2000);
data2.var = randn(63,2000);
data2.dof = randn(63,2000);
data2.time = (1:2000)/1000;
for i=1:63
  data2.label{i} = num2str(i);
end
data2.trial = randn(42,63,2000);
data2.trialinfo = randn(42,5);
data2.dimord = 'rpt_chan_time';

% the following is known to fail with the latest version on 8 April 2014
% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=2528

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
output = ft_math(cfg, data1, data2);

