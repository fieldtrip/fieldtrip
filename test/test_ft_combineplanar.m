function test_ft_combineplanar

% WALLTIME 00:20:00
% MEM 1gb
% DEPENDENCY ft_combineplanar ft_preprocessing ft_timelockanalysis ft_prepare_neighbours ft_megplanar
% DATA no

headmodel = [];
headmodel.r = 0.09;
headmodel.o = [0 0 0.04];
headmodel.unit = 'm';

%%

% these are in the fieldtrip/template/gradiometer directory
senstype = {
  'ctf64.mat'
  'ctf151.mat'
  'ctf275.mat'
  'bti148.mat'
  'bti248.mat'
  'itab153.mat'
  'yokogawa160.mat'
  };

for i=1:numel(senstype)
  grad = ft_read_sens(senstype{i});

  cfg = [];
  cfg.headmodel = headmodel;
  cfg.grad = grad;
  data = ft_dipolesimulation(cfg);

  cfg = [];
  timelock = ft_timelockanalysis(cfg, data);

  cfg              = [];
  cfg.feedback     = 'no';
  cfg.method       = 'template';
  cfg.neighbours   = ft_prepare_neighbours(cfg, timelock);
  cfg.planarmethod = 'sincos';
  timelock_planar = ft_megplanar(cfg, timelock);

  cfg = [];
  cfg.updatesens = 'no';
  timelock_combined = ft_combineplanar(cfg, timelock_planar);

end

%%

% these are in the fieldtrip/template/gradiometer directory
senstype = {
  'neuromag122.mat'
  'neuromag306.mat'
  };

for i=1:numel(senstype)
  grad = ft_read_sens(senstype{i});

  cfg = [];
  cfg.headmodel = headmodel;
  cfg.grad = grad;
  data_planar = ft_dipolesimulation(cfg);

  cfg = [];
  timelock_planar = ft_timelockanalysis(cfg, data_planar);

  cfg = [];
  cfg.updatesens = 'no';
  timelock_combined = ft_combineplanar(cfg, timelock_planar);

end