function test_pull929

% MEM 1000mb
% WALLTIME 00:10:00

% TEST ft_scalpcurrentdensity ft_fetch_sens

%% load data
% the data contains an ERP with 128 channels and 2001 time bins
% the ERP has a fronto-central negativity at 100ms
% the SCD resolves into two bilateral blobs

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2685/bug2685.mat'));

ERP_standard = rmfield(ERP_standard, 'grad');

cfg = [];
cfg.center = 'yes';
layout = ft_prepare_layout(cfg, ERP_standard);

cfg = [];
cfg.layout = layout;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, ERP_standard);


%% 1a) Test method spline in absence of bad channels

gooddata = ERP_standard;
cfg = [];
cfg.method = 'spline';
goodscd = ft_scalpcurrentdensity(cfg, gooddata);

cfg = [];
cfg.layout = layout;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, goodscd);


%% 1b) Test method spline in presence of bad channels

baddata = ERP_standard;
baddata.avg(1:64,:) = NaN; % set first half of channels to NaN

cfg = [];
cfg.method = 'spline';
badscd = ft_scalpcurrentdensity(cfg, baddata);

cfg = [];
cfg.layout = layout;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, goodscd, badscd);

% Compare difference between scd with and without bad channels:
diffscd = goodscd.avg - badscd.avg;
if mean(abs(diffscd(:))) > 1e-08 % even with 50% of channels missing, this threshold should not be crossed
  error('scd with and without bad channels strongly differs');
end


%% 1c) Test method spline with constant offset

offsetdata = ERP_standard;
offsetdata.avg = gooddata.avg+1000; % add constant offset of 1000 V
cfg = [];
cfg.method = 'spline';
offsetscd = ft_scalpcurrentdensity(cfg, offsetdata);

cfg = [];
cfg.layout = layout;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, goodscd, offsetscd);

% Compare difference between scd with and without offset:
diffscd = goodscd.avg - offsetscd.avg;
if mean(abs(diffscd(:))) > 100*eps % only floating point errors expected
  error('scd with and without offset strongly differs');
end


%% 1d) Test method spline with data channels shuffled

shuffledata = ERP_standard;
randidx = randperm(length(shuffledata.label));
[~,unidx] = sort(randidx); % unshuffle
shuffledata.avg = shuffledata.avg(randidx,:);
shuffledata.var = shuffledata.var(randidx,:);
shuffledata.dof = shuffledata.dof(randidx,:);
shuffledata.label = shuffledata.label(randidx);

cfg = [];
cfg.method = 'spline';
shufflescd = ft_scalpcurrentdensity(cfg, shuffledata);

cfg = [];
cfg.layout = layout;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, goodscd, shufflescd);

% Compare difference between scd with and without channels shuffled:
diffscd = goodscd.avg - shufflescd.avg(unidx,:);
if mean(abs(diffscd(:))) > 100*eps % only floating point errors expected
  error('scd with and without shuffling of channels strongly differs');
end


%% 1e) Test method spline with electrodes shuffled

shuffledata = ERP_standard;
shuffledata.elec.chanpos = shuffledata.elec.chanpos(randidx,:);
shuffledata.elec.elecpos = shuffledata.elec.elecpos(randidx,:);
shuffledata.elec.label = shuffledata.elec.label(randidx);

cfg = [];
cfg.method = 'spline';
shufflescd = ft_scalpcurrentdensity(cfg, shuffledata);

cfg = [];
cfg.layout = layout;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, goodscd, shufflescd);

% Compare difference between scd with and without channels shuffled:
diffscd = goodscd.avg - shufflescd.avg;
if mean(abs(diffscd(:))) > 100*eps % only floating point errors expected
  error('scd with and without shuffling of channels strongly differs');
end


%% 2a) Test method finite in absence of bad channels

gooddata = ERP_standard;

cfg = [];
cfg.method = 'finite';
scd = ft_scalpcurrentdensity(cfg, gooddata);


%% 2b) Test method finite in presence of bad channels

baddata = ERP_standard;
baddata.avg(1:64,:) = NaN; % set first half of channels to NaN

try
  cfg = [];
  cfg.method = 'finite';
  failed = true;
  scd = ft_scalpcurrentdensity(cfg, baddata); % should fail on data with channels containing NaNs
  failed = false;
end

if failed
  fprintf('this correctly throws error in presence of bad channels\n')
else
  error('this should have thrown an error in presence of bad channels\n')
end
