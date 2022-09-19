function test_ft_connectivitysimulation

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_connectivitysimulation

cfg = [];
cfg.method      = 'ar';
cfg.nsignal     = 32;
cfg.ntrials     = 100;
cfg.triallength = 2;
cfg.fsample     = 500;

% test 'ar' method
cfg.params = zeros(cfg.nsignal,cfg.nsignal,cfg.ntrials);
cfg.params(14,21,42) = 0.5;
cfg.noisecov = 0.5*eye(cfg.nsignal);
dataout = ft_connectivitysimulation(cfg);

return
%% test linear_mix method => wait for issue to be solved
cfg.mix = 0.5*ones(32,2);
cfg.delay = 50*ones(32,2);

dataout = ft_connectivitysimulation(cfg);
