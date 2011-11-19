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
data2 = ft_redefinetrial(cfg, data);
isfield(data2, 'cfg')

cfg = [];
cfg.begsample = 10;
cfg.endsample = 20;
data2 = ft_redefinetrial(cfg, data);
isfield(data2, 'cfg')