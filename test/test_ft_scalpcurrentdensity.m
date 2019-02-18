function test_ft_scalpcurrentdensity

% MEM 1000mb
% WALLTIME 00:01:00

% TEST ft_scalpcurrentdensity ft_fetch_sens

%% load data
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2685/bug2685.mat'));
% 128 channels, 1 trial with 2001 time bins
%% Plot ERP of unfiltered data
figure;
plot(ERP_standard.avg')
title('ERP on unfiltered data');
axis([0 size(ERP_standard.avg,2) min(ERP_standard.avg(:)) max(ERP_standard.avg(:))])
%% 1a) Test method spline in absence of bad channels
gooddata            = ERP_standard;
cfg                 = [];
cfg.method          = 'spline';
cfg.elec            = gooddata.elec;
goodscd   = ft_scalpcurrentdensity(cfg, gooddata);
% Plot scd:
figure;
plot(goodscd.avg')
title('ERP on scd without bad channels');
axis([0 size(goodscd.avg,2) min(goodscd.avg(:)) max(goodscd.avg(:))])
%% 1b) Test method spline in presence of bad channels
baddata             = ERP_standard;
baddata.avg(1:2,1) = NaN; % set first sample of first half of channels to NaN
cfg                 = [];
cfg.method          = 'spline';
cfg.elec            = baddata.elec;
badscd   = ft_scalpcurrentdensity(cfg, baddata);
% Plot scd:
figure;
plot(badscd.avg')
title('ERP on scd with bad channels');
%axis([0 size(badscd.avg,2) min(badscd.avg(:)) max(badscd.avg(:))])
axis([0 size(goodscd.avg,2) min(goodscd.avg(:)) max(goodscd.avg(:))])
% Compare difference between scd with and without bad channels:
diffscd = goodscd.avg - badscd.avg;
if mean(abs(diffscd(:))) > 1e-08 % even with 50% of channels missing, this threshold should not be crossed
    error('scd with and without bad channels strongly differs');
end
%% 1c) Test method spline with constant offset
offsetdata          = ERP_standard;
offsetdata.avg      = gooddata.avg+1000; % add constant offset of 1000 V
cfg                 = [];
cfg.method          = 'spline';
cfg.elec            = offsetdata.elec;
offsetscd   = ft_scalpcurrentdensity(cfg, offsetdata);
% Plot scd:
figure;
plot(offsetscd.avg')
title('ERP on scd with offset');
%axis([0 size(offsetscd.avg,2) min(offsetscd.avg(:)) max(offsetscd.avg(:))])
axis([0 size(goodscd.avg,2) min(goodscd.avg(:)) max(goodscd.avg(:))])
% Compare difference between scd with and without offset:
diffscd = goodscd.avg - offsetscd.avg;
if mean(abs(diffscd(:))) > 1e-15
    error('scd with and without offset strongly differs');
end
%% 1d) Test method spline with channels shuffled
shuffledata         = ERP_standard;
randidx             = randperm(length(shuffledata.label));
[~,unidx]           = sort(randidx); % unshuffle
shuffledata.avg     = shuffledata.avg(randidx,:);
shuffledata.var     = shuffledata.var(randidx,:);
shuffledata.dof     = shuffledata.dof(randidx,:);
shuffledata.label   = shuffledata.label(randidx);
shuffledata.elec.chanpos   = shuffledata.elec.chanpos(randidx,:);
shuffledata.elec.elecpos   = shuffledata.elec.elecpos(randidx,:);
shuffledata.elec.label     = shuffledata.elec.label(randidx);
shuffledata.cfg.channel    = shuffledata.cfg.channel(randidx);
cfg                 = [];
cfg.method          = 'spline';
cfg.elec            = shuffledata.elec;
shufflescd   = ft_scalpcurrentdensity(cfg, shuffledata);
% Plot scd:
figure;
plot(shufflescd.avg(unidx,:)')
title('ERP on scd with channels shuffled');
%axis([0 size(shuffledata.avg,2) min(shuffledata.avg(:)) max(shuffledata.avg(:))])
axis([0 size(goodscd.avg,2) min(goodscd.avg(:)) max(goodscd.avg(:))])
% Compare difference between scd with and without channels shuffled:
diffscd = goodscd.avg - shufflescd.avg(unidx,:);
if mean(abs(diffscd(:))) > 1e-19 % only floating point errors expected
    error('scd with and without shuffling of channels strongly differs');
end
%% 2a) Test method finite in absence of bad channels
gooddata            = ERP_standard;
cfg                 = [];
cfg.method          = 'finite';
cfg.elec            = ERP_standard.elec;
scd   = ft_scalpcurrentdensity(cfg, gooddata);
%% 2b) Test method finite in presence of bad channels
try
    failed = true;
    baddata             = ERP_standard;
    cfg                 = [];
    cfg.method          = 'finite';
    cfg.elec            = ERP_standard.elec;
    scd   = ft_scalpcurrentdensity(cfg, baddata); % should fail on data with channels containing NaNs
    failed = false;
end
if failed
    fprintf('method finite correctly throws error in presence of bad channels\n')
else
    error('Method finite should throw error if some channels bad, but does not\n')
end
