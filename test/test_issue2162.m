function test_issue2162

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_checkconfig ft_postamble_history
% DATA no

%%

nchan = 10;
ntrial = 5;
ntime = 1000;
fsample = 1000;

data = [];
for i=1:ntrial
  data.time{i} = (0:ntime-1)/fsample;
  data.trial{i} = randn(nchan, ntime);
end
for i=1:nchan
  data.label{i} = sprintf('%d', i);
end

%%

cfg = [];
cfg.something = zeros(1e3,1e3);
cfg.checksize = 1e6;

cfg = ft_checkconfig(cfg, 'checksize', 'yes');

assert(ischar(cfg.something)) % should be cleared by checkconfig

%%

cfg = [];
cfg.something = zeros(1e3,1e3);
cfg.checksize = 1e6;

timelock = ft_timelockanalysis(cfg, data);

assert(ischar(timelock.cfg.something)) % should be cleared by checkconfig

