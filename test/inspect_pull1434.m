function inspect_pull1434

% WALLTIME 00:10:00
% MEM 1gb
% DATA no

%%

fsample = 1000;
nsamples = 10*fsample;
nchans = 10;

data = [];
for i=1:nchans
  data.label{i} = num2str(i);
end
data.trial{1} = randn(nchans, nsamples);
data.time{1} = (1:nsamples)/fsample;


nevent = 300;
event = [];
for i=1:nevent
  event(i).type = 'Trigger';
  event(i).value = i;
  event(i).sample = (i-1) * round(nsamples/nevent) + 1;
end

%%


cfg = [];
cfg.event = event;
cfg.viewmode = 'vertical';
cfg.ylim = [-5 5];
ft_databrowser(cfg, data);
