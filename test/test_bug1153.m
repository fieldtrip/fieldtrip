function test_bug1153

% TEST test_bug1153 ft_redefinetrial

% problem: ft_redefinetrial loses the cfg in the output when specifying
% cfg.trl

data = [];
data.trial{1} = randn(1,100);
data.time{1}  = 1:100;
data.label{1} = 'chan1';
data.sampleinfo = [1 100];

cfg = [];
cfg.trl = [10 20 0];
data = ft_redefinetrial(cfg, data);
isfield(data, 'cfg')