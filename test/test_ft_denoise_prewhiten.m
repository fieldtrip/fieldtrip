function test_ft_denoise_prewhiten

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_denoise_prewhiten

% create some data
data = [];
data.label = {'chan01';'chan02';'chan03'};

% full rank data
mix = [1 0 0;0 1 0;0.5 0.5 0.1];
for k = 1:10
  data.trial{k} = mix*randn(3,1000);
  data.time{k}  = (0:999)./1000;
end

cfg = [];
cfg.covariance = 'yes';
tlck = ft_timelockanalysis(cfg, data);

cfg = [];
datawhite = ft_denoise_prewhiten(cfg, data, tlck);

cfg = [];
cfg.covariance = 'yes';
tlckwhite = ft_timelockanalysis(cfg, datawhite);

cfg = [];
tlckwhite2 = ft_denoise_prewhiten(cfg, tlck, tlck);

% rank deficient data
mix = [1 0 0;0 1 0;0.5 0.5 0];
for k = 1:10
  data.trial{k} = mix*randn(3,1000);
  data.time{k}  = (0:999)./1000;
end

cfg = [];
cfg.covariance = 'yes';
tlck = ft_timelockanalysis(cfg, data);

cfg = [];
datawhite = ft_denoise_prewhiten(cfg, data, tlck);

cfg = [];
cfg.covariance = 'yes';
tlckwhite = ft_timelockanalysis(cfg, datawhite);

cfg = [];
tlckwhite2 = ft_denoise_prewhiten(cfg, tlck, tlck);

