function test_bug2647

% TEST ft_freqstatistics

% WALLTIME 00:10:00
% MEM 1gb

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2647.mat'));
stat = ft_freqstatistics(cfg,data,BL);
%assert(all(~isfinite(stat.stat(:))));
assert(~all(~isfinite(stat.stat(:)))); % this should now work after making ft_statfun_actvsblT more robust for NaNs

cfg.design = cfg.design(:,[1:4 6:9]);
data.powspctrm = data.powspctrm(2:end,:,:,:);
BL.powspctrm   = BL.powspctrm(2:end,:,:,:);
stat = ft_freqstatistics(cfg,data,BL);
assert(~all(~isfinite(stat.stat(:))));

% conclusion: the first rpt was all NaN, causing stat.stat to be all NaN,
% no meaningful inference possible
