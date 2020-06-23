%function test_ft_eventtiminganalysis

% MEM 2gb
% WALLTIME 00:15:00
% DEPENDENCY ft_eventtiminganalysis

load(dccnpath('/home/common/matlab/fieldtrip/data/test/test_ft_eventtiminganalysis.mat'));
ft_hastoolbox('lagextraction', 1);

cfg = [];
cfg.method = 'gbve';
cfg.gbve.latency = [0 inf];
cfg.gbve.use_maximum = false;
out_gbve = ft_eventtiminganalysis(cfg, data);

cfg = [];
cfg.method = 'aseo';
cfg.aseo.initlatency = [0.2 0.5];
cfg.aseo.jitter = 0.15;
out_aseo = ft_eventtiminganalysis(cfg, data);

rt=(data.trialinfo(:,end)-data.trialinfo(:,end-1))./1200;
figure;plot(rt, out_gbve.params.latency,'o');
hold on;plot(rt, out_aseo.params.latency+mean(rt),'o');
