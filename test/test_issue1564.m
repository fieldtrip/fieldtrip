function test_issue1564

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_artifact_clip ft_artifact_threshold ft_artifact_jump ft_rejectartifact

global ft_default
ft_default.representation = 'table';

%%

ntrials = 10;
nchans = 3;
fsample = 1000;

data = [];
data.label = arrayfun(@num2str, 1:nchans, 'UniformOutput', false);
for i=1:ntrials
  data.time{i} = (1:fsample)/fsample;
  data.trial{i} = randn(nchans, fsample);
  data.sampleinfo(i,1) = (i-1)*fsample + 1 + (i*100); % create some gaps between the trials
  data.sampleinfo(i,2) = data.sampleinfo(i,1) + fsample - 1;
end

if strcmp(ft_default.representation, 'numeric')
  data.trialinfo = [1:ntrials]';
else
  trialnum = arrayfun(@num2str, 1:ntrials, 'UniformOutput', false)';
  data.trialinfo = table(trialnum);
end

%%

cfg = [];
cfg.continuous = 'no';
ft_databrowser(cfg, data);

%%

cfg = [];
cfg.continuous = 'yes';
ft_databrowser(cfg, data);

%%

cfg = [];
cfg.continuous = [];
ft_databrowser(cfg, data);


%%
% close all

data_clip = data;
data_clip.trial{1}(1,201:700) =  2;
data_clip.trial{2}(2,301:400) = -2;
data_clip.trial{2}(2,501:600) = -2;

cfg = [];
cfg.artfctdef.clip.amplthreshold = 0;
[cfg, artifact] = ft_artifact_clip(cfg, data_clip);

% cfg = [];
% cfg.artfctdef.clip.amplthreshold = '0%';
% [cfg, artifact] = ft_artifact_clip(cfg, data_clip)

% cfg = [];
% cfg.artfctdef.clip.amplthreshold = '100%';
% [cfg, artifact] = ft_artifact_clip(cfg, data_clip)

% explore the detected artifacts in the databrowser
ft_databrowser(cfg, data_clip);

%%

cfg.artfctdef.reject = 'partial';
cfg.artfctdef.minaccepttim = 0;
cfg.artfctdef.feedback = 'yes';
data_rej = ft_rejectartifact(cfg, data_clip);
if ~istable(data.trialinfo)
  assert(isequal(data_rej.trialinfo, [1 1 2 2 2 3 4 5 6 7 8 9 10]'));
end

cfg.artfctdef.reject = 'complete';
data_rej = ft_rejectartifact(cfg, data_clip);
if ~istable(data.trialinfo)
  assert(isequal(data_rej.trialinfo, [3 4 5 6 7 8 9 10]'));
end

cfg.artfctdef.reject = 'nan';
data_rej = ft_rejectartifact(cfg, data_clip);
if ~istable(data.trialinfo)
  assert(isequal(data_rej.trialinfo, [1 2 3 4 5 6 7 8 9 10]'));
end


%%
close all

data_thresh = data;
data_thresh.trial{1}(1,201:700) = data_thresh.trial{1}(1,201:700) + 11;
data_thresh.trial{2}(2,301:400) = data_thresh.trial{2}(2,301:400) + 12;
data_thresh.trial{2}(3,501:600) = data_thresh.trial{2}(3,501:600) - 13;

cfg = [];
cfg.artfctdef.threshold.channel = 'all';
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.min = -8;
cfg.artfctdef.threshold.max =  8;

[cfg, artifact] = ft_artifact_threshold(cfg, data_thresh);

% explore the detected artifacts in the databrowser
ft_databrowser(cfg, data_thresh);

cfg.artfctdef.reject = 'partial';
cfg.artfctdef.minaccepttim = 0;
cfg.artfctdef.feedback = 'yes';
data_rej = ft_rejectartifact(cfg, data_thresh);
if ~istable(data.trialinfo)
  assert(isequal(data_rej.trialinfo, [1 1 2 2 2 3 4 5 6 7 8 9 10]'));
end

cfg.artfctdef.reject = 'complete';
data_rej = ft_rejectartifact(cfg, data_thresh);
if ~istable(data.trialinfo)
  assert(isequal(data_rej.trialinfo, [3 4 5 6 7 8 9 10]'));
end

cfg.artfctdef.reject = 'nan';
data_rej = ft_rejectartifact(cfg, data_thresh);
if ~istable(data.trialinfo)
  assert(isequal(data_rej.trialinfo, [1 2 3 4 5 6 7 8 9 10]'));
end


%%
% Although I identified in issue 1564 that the detected jump artifacts could also
% have their channel in the outpout table, in the end I did not implement it. It
% would be a lot of work and not really worth it.

close all

data_jump = data;
data_jump.trial{1}(1,401:end) = data_jump.trial{1}(1,401:end) + 15;
data_jump.trial{2}(2,601:end) = data_jump.trial{2}(2,601:end) - 15;

cfg = [];
cfg.artfctdef.jump.channel = 'all';
cfg.artfctdef.jump.trlpadding = 0;
cfg.artfctdef.jump.fltpadding = 0;
cfg.artfctdef.jump.artpadding = 0.1;
cfg.artfctdef.jump.cutoff = 20;
cfg.artfctdef.jump.interactive = 'no';

[cfg, artifact] = ft_artifact_jump(cfg, data_jump);

% explore the detected artifacts in the databrowser
ft_databrowser(cfg, data_jump);
