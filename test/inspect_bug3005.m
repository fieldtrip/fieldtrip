function inspect_bug3005

%%

nchan = 10;
ntime = 1000;
fsample = 1000;
ntrial = 1; % this is important

data = [];
for i=1:nchan
  data.label{i} = sprintf('%d', i);
end
for i=1:ntrial
  data.trial{i} = randn(nchan, ntime);
  data.time{i} = (1:ntime)/fsample;
end

%%

cfg = [];
cfg.method = 'summary';
rej1 = ft_rejectvisual(cfg, data);

%%

cfg = [];
cfg.method = 'trial';
rej2 = ft_rejectvisual(cfg, data);

%%

cfg = [];
cfg.method = 'channel';
rej3 = ft_rejectvisual(cfg, data);

