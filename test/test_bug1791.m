function test_bug1791

% TEST test_bug1791
% TEST ft_removetmsartifact ft_artifact_tmspulse

load /home/common/matlab/fieldtrip/data/test/bug1791/data.mat;

cfg = [];
cfg.pulsewidth = 0.010;
cfg.offset = -0.005;
cfg.method = 'twopassfilter';
cfg.latency = 0.008/2;
data1 = ft_removetmsartifact(cfg, data);

cfg = [];
cfg.method = 'interpolatepulse';
cfg.latency = 0.008/2;
cfg.pulsewidth = 0.010;
cfg.offset = -0.005;
data2 = ft_removetmsartifact(cfg, data);

figure 
hold on
plot(data.time{2},  data.trial{2}(1,:),  'b')
plot(data1.time{2}, data1.trial{2}(1,:), 'r')
plot(data2.time{2}, data2.trial{2}(1,:), 'g')
legend({'original';'twopass';'interpolate'});

