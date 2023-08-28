function test_bug2783

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_redefinetrial ft_checkdata
% DATA no

%% create some data

data = [];
data.time = 1:1000;
data.trial = ones(1,3,1000);
data.dimord = 'rpt_chan_time';
data.label = {'a','b','c'};

%% redefine

cfg = [];
cfg.offset = 100;
dat2 = ft_redefinetrial(cfg, data);

%% check

assert(isfield(dat2, 'trial') && isequal(dat2.trial, data.trial));

end
