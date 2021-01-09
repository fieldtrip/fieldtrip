function test_bug3238

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

% this test script is under development and is known to fail at the moment
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3238

% test the handling of partial artifacts with nans

%% create some raw data

nchan = 16;
ntrial = 5;
ntime = 10000;
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
data.trial{2}(:,1001:1100) = nan;

% 3nd trial, short time segment bad in one channel
data.trial{3}(2,1201:1300) = nan;

% 4th trial, completely bad;
data.trial{4}(:,:) = nan;

%% check the nan-awareness of the various ft_preproc_ functions
label = cell(1,0);
for k = 1:ntrial
  datout = ft_preproc_bandpassfilter(data.trial{k}, fsample, [15 25], [], 'firws');
  nnan(:,1,k) = sum(isnan(datout),2);
  label{end+1} = 'bandpassfilter';
  datout = ft_preproc_bandstopfilter(data.trial{k}, fsample, [15 25], [], 'firws');
  nnan(:,2,k) = sum(isnan(datout),2);
  label{end+1} = 'bandstopfilter';
  datout = ft_preproc_baselinecorrect(data.trial{k});
  nnan(:,3,k) = sum(isnan(datout),2);
  label{end+1} = 'baselinecorrect';
  datout = ft_preproc_denoise(data.trial{k}, data.trial{5});
  nnan(:,4,k) = sum(isnan(datout),2);
  label{end+1} = 'denoise';
  datout = ft_preproc_derivative(data.trial{k});
  nnan(:,5,k) = sum(isnan(datout),2);
  label{end+1} = 'derivative';
  datout = ft_preproc_detrend(data.trial{k});
  nnan(:,6,k) = sum(isnan(datout),2);
  label{end+1} = 'detrend';
  datout = ft_preproc_dftfilter(data.trial{k}, fsample, 50);
  nnan(:,7,k) = sum(isnan(datout),2);
  label{end+1} = 'dftfilter';
  datout = ft_preproc_highpassfilter(data.trial{k}, fsample, 10, [], 'firws');
  nnan(:,8,k) = sum(isnan(datout),2);
  label{end+1} = 'highpassfilter';
  datout = ft_preproc_hilbert(data.trial{k});
  nnan(:,9,k) = sum(isnan(datout),2);
  label{end+1} = 'hilbert';
  datout = ft_preproc_lowpassfilter(data.trial{k}, fsample, 40, [], 'firws');
  nnan(:,10,k) = sum(isnan(datout),2);
  label{end+1} = 'lowpassfilter';
  datout = ft_preproc_medianfilter(data.trial{k}, 5);
  nnan(:,11,k) = sum(isnan(datout),2);
  label{end+1} = 'medianfilter';
  %ft_preproc_online_downsample_apply
  %ft_preproc_online_downsample_init
  %ft_preproc_online_filter_apply
  %ft_preproc_online_filter_init
  datout = ft_preproc_padding(data.trial{k}, 'mean', 50);
  nnan(:,12,k) = sum(isnan(datout),2);
  label{end+1} = 'padding';
  datout = ft_preproc_polyremoval(data.trial{k}, 2);
  nnan(:,13,k) = sum(isnan(datout),2);
  label{end+1} = 'polyremoval';
  datout = ft_preproc_rectify(data.trial{k});
  nnan(:,14,k) = sum(isnan(datout),2);
  label{end+1} = 'rectify';
  datout = ft_preproc_rereference(data.trial{k}, 1:nchan, 'avg', true);
  nnan(:,15,k) = sum(isnan(datout),2);
  label{end+1} = 'rereference';
  datout = ft_preproc_resample(data.trial{k}, fsample, 500, 'resample');
  nnan(:,16,k) = sum(isnan(datout),2);
  label{end+1} = 'resample';
  datout = ft_preproc_slidingrange(data.trial{k}, 5);
  nnan(:,17,k) = sum(isnan(datout),2);
  label{end+1} = 'slidingrange';
  datout = ft_preproc_smooth(data.trial{k}, 5);
  nnan(:,18,k) = sum(isnan(datout),2);
  label{end+1} = 'smooth';
  datout = ft_preproc_standardize(data.trial{k});
  nnan(:,19,k) = sum(isnan(datout),2);
  label{end+1} = 'standardize';
end

figure; imagesc(nnan(:,:,1)');
set(gca, 'ytick', 1:19, 'yticklabel', label); title('trial 1'); xlabel('channel');

figure; imagesc(nnan(:,:,2)');
set(gca, 'ytick', 1:19, 'yticklabel', label); title('trial 2'); xlabel('channel');

figure; imagesc(nnan(:,:,3)');
set(gca, 'ytick', 1:19, 'yticklabel', label); title('trial 3'); xlabel('channel');

figure; imagesc(nnan(:,:,4)');
set(gca, 'ytick', 1:19, 'yticklabel', label); title('trial 4'); xlabel('channel');

figure; imagesc(nnan(:,:,5)');
set(gca, 'ytick', 1:19, 'yticklabel', label); title('trial 5'); xlabel('channel');

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

cfg = [];
cfg.keeptrials = 'no';
timelock2b = ft_timelockanalysis(cfg, data);
assert(~any(isnan(timelock2b.avg(:)))); % but not here
assert(~any(isnan(timelock2b.var(:)))); % or here

% TODO check the df or the dof field (if present)

%%

% cfg = [];
% cfg.covariance = 'yes';
% timelock3 = ft_timelockanalysis(cfg, data);
% 
% assert(~any(isnan(timelock3.cov(:))));

% TODO check the df or the dof field (if present)

%%

% cfg = [];
% cfg.keeptrials = 'yes';
% cfg.covariance = 'yes';
% timelock4 = ft_timelockanalysis(cfg, data);
% 
% tmp = squeeze(timelock4.cov(1,:,:)); assert(~any(isnan(tmp(:))));
% tmp = squeeze(timelock4.cov(2,:,:)); assert(~any(isnan(tmp(:))));
% tmp = squeeze(timelock4.cov(3,:,:)); assert(~any(isnan(tmp(:))));
% tmp = squeeze(timelock4.cov(4,:,:)); assert( all(isnan(tmp(:)))); % this is all nan

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



