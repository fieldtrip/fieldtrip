function failed_bug3238

% WALLTIME 00:10:00
% MEM 2gb

% this test script is under development and is known to fail at the moment
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3238

% test the handling of partial artifacts with nans

%% create some raw data

nchan = 16;
ntrial = 10;
ntime = 1000;
fsample = 1000;

data = [];
for i=1:nchan
  data.label{i,1} = num2str(i);
end
for i=1:ntrial
  data.trial{i} = randn(nchan,ntime);
  data.time{i} = (1:ntime)/fsample;
end

% 1st trial, channel 1 is bad
data.trial{1}(1,:) = nan;

% 2nd trial, time segment is bad for all channels
data.trial{2}(:,1:100) = nan;

% 3nd trial, short time segment bad in one channel
data.trial{3}(2,201:300) = nan;

% 4th trial, completely bad;
data.trial{4}(:,:) = nan;

%%

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 50];
data_flt1 = ft_preprocessing(cfg, data);

assert(sum(any(isnan(data_flt1.trial{1}), 2))==1); % one channel
assert(sum(any(isnan(data_flt1.trial{3}), 2))==1); % one channel

assert(sum(all(isnan(data_flt1.trial{2}), 2))==nchan); % all channels
assert(sum(all(isnan(data_flt1.trial{4}), 2))==nchan); % all channels

%%

cfg = [];
cfg.reref = 'yes';
cfg.refchannel = 'all';
data_flt2 = ft_preprocessing(cfg, data);

assert(sum(any(isnan(data_flt2.trial{1}), 2))==nchan); % one channel
assert(sum(any(isnan(data_flt2.trial{3}), 2))==nchan); % one channel

assert(sum(any(isnan(data_flt2.trial{2}), 2))==nchan); % all channels
assert(sum(any(isnan(data_flt2.trial{4}), 2))==nchan); % all channels

%%

cfg = [];
cfg.keeptrials = 'no';
timelock1 = ft_timelockanalysis(cfg, data);

assert(~any(isnan(timelock1.avg(:))));
assert(~any(isnan(timelock1.var(:))));

% TODO check the df or the dof field (if present)

%%

cfg = [];
cfg.keeptrials = 'yes';
timelock2 = ft_timelockanalysis(cfg, data);

assert(any(isnan(timelock2.trial(:)))); % there should be nans in here

assert(~any(isnan(timelock2.avg(:)))); % but not here
assert(~any(isnan(timelock2.var(:)))); % or here

% TODO check the df or the dof field (if present)

%%

cfg = [];
cfg.covariance = 'yes';
timelock3 = ft_timelockanalysis(cfg, data);

assert(~any(isnan(timelock3.cov(:))));

% TODO check the df or the dof field (if present)

%%

cfg = [];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes';
timelock4 = ft_timelockanalysis(cfg, data);

tmp = squeeze(timelock4.cov(1,:,:)); assert(~any(isnan(tmp(:))));
tmp = squeeze(timelock4.cov(2,:,:)); assert(~any(isnan(tmp(:))));
tmp = squeeze(timelock4.cov(3,:,:)); assert(~any(isnan(tmp(:))));
tmp = squeeze(timelock4.cov(4,:,:)); assert( all(isnan(tmp(:)))); % this is all nan

% TODO check the df or the dof field (if present)


%%

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'hanning';
cfg.keeptrials = 'yes';
freq1 = ft_freqanalysis(cfg, data);

% todo: look at output powspctrm

%%

cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.taper = 'hanning';
cfg.foi = 1:30;
cfg.toi = 0:0.05:1;
cfg.t_ftimwin = 0.5 * ones(size(cfg.foi));
% cfg.keeptrials = 'yes';
freq2 = ft_freqanalysis(cfg, data);

% todo: look at output powspctrm, should be nan, regardless of time

%%

cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.foilim = [1 30];
cfg.toi = 0:0.05:1;
% cfg.keeptrials = 'yes';
freq3 = ft_freqanalysis(cfg, data);

% todo: look at output powspctrm

%%

% todo
% - ft_multiplotER
% - ft_topoplotER
% - ft_rejectvisual
% - ft_databrowser

