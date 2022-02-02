function inspect_pull1946

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_databrowser

%%
% create some uniform data

nchan = 10;
ntrial = 30;
nsample = 1000;
fsample = 1000;

data = [];
for i=1:nchan
  data.label{i} = num2str(i);
end

for i=1:ntrial
  data.trial{i} = randn(nchan,nsample);
  data.time{i}  = (1:nsample)/fsample;
end

data.sampleinfo(:,1) = ((1:ntrial)-1)*nsample+1;
data.sampleinfo(:,2) = ((1:ntrial)  )*nsample;

% add an artifact to channel 2, trial 2
data.trial{2}(2,:) = 3 * data.trial{2}(2,:);

%%

cfg = [];
cfg.method = 'summary';
cfg.keeptrials = 'yes';
cfg.keepchannels = 'yes';
cfg.channel = setdiff(1:nchan, 2);
cfg.trials = setdiff(1:ntrial, 2);
dataclean = ft_rejectvisual(cfg, data);


