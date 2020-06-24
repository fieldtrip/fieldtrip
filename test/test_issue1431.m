function test_issue1431

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_heartrate ft_interpolatenan

%%

qrs = gradient(gradient(chebwin(102)));
qrs = qrs(2:end-1);
qrs = -qrs + qrs(1);
qrs = qrs/max(qrs);
% plot(qrs);

%%

nchan = 1;
ntrl = 60;
ntime = 1000;
fsample = 1000;

data = [];
data.label = {'ECG'};
for i=1:ntrl
  % here we put a beat in the middle of each 1-second segment, but with some variation
  
  dat = zeros(nchan,ntime);
  sel = (451:550) + round(50*randn(1));
  dat(sel) = qrs;
  
  data.time{i} = ((i-1)*ntime + (1:ntime))/fsample;
  data.trial{i} = dat + 0.01*randn(size(dat));
end

% concatenate all segments to get a single continuous representation with 60 bpm
data.trial = {cat(2, data.trial{:})};
data.time  = {cat(2, data.time{:})};

%%

cfg = [];
cfg.blocksize = 60;
ft_databrowser(cfg, data);

%%

cfg = [];
cfg.threshold = 2;
cfg.method = 'findpeaks';
cfg.channel = 'ECG';
heartrate = ft_heartrate(cfg, data);

% there should be 60 beats
assert(sum(heartrate.trial{1}(strcmp(heartrate.label, 'heartbeatonset'),:))==60);

cfg = [];
cfg.method = 'pantompkin';
cfg.channel = 'ECG';
heartrate = ft_heartrate(cfg, data);

% there should be 60 beats
assert(sum(heartrate.trial{1}(strcmp(heartrate.label, 'heartbeatonset'),:))==60);

%%

cfg = [];
% find the nans at the edges
cfg = ft_artifact_nan(cfg, heartrate);

% insert some manual artifacts
cfg.artfctdef.manual.artifact = [
  5001  15000
  25001 26000
  ];

cfg.blocksize = 60;
cfg.channel = {'heartrate'};
ft_databrowser(cfg, heartrate);

% this replaces the manually flagged pieces of data with nan
cfg.reject = 'nan';
heartrate_rej = ft_rejectartifact(cfg, heartrate);

figure
plot(heartrate_rej.time{1}, heartrate_rej.trial{1}(1,:))

%%

close all

method = {
  'linear'
  'nearest'
  'spline'
  'pchip'
  'cubic'
  'v5cubic'
  'makima'
  };

figure
plot(heartrate_rej.time{1}, heartrate_rej.trial{1}(1,:))
title('original');

for i=1:numel(method)
  cfg = [];
  cfg.prewindow = 0.005;
  cfg.postwindow = 0.005;
  cfg.method = method{i};
  heartrate_interp = ft_interpolatenan(cfg, heartrate_rej);
  
  % there should not be any nan left
  assert(~any(isnan(heartrate_interp.trial{1}(:))));
  
  figure
  plot(heartrate_interp.time{1}, heartrate_interp.trial{1}(1,:))
  title(method{i});
end
