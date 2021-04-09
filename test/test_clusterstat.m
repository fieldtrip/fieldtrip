function test_clusterstat

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY findcluster clusterstat ft_statistics_montecarlo

%%%%%%%%5
% some code -> this relates to issue 1732
pattern = [...
    1 1 1 1 0 1 1 1 1 1;...
    1 0 0 0 0 0 0 1 0 0;...
    1 1 1 0 0 0 0 1 0 0;...
    1 0 0 0 0 0 0 1 0 0;...
    1 0 0 0 0 0 0 1 0 0;...
    ]; % FT

%% 0D (per channel)
% make fake dataset
data = cell(1,10);
for idat = 1:10
  data{idat}.label = {'ch1' 'ch2' 'ch3' 'ch4' 'ch5'}';
  data{idat}.dimord = 'chan_freq';
  data{idat}.parameter = rand(5,1);
  if idat <= 5
    data{idat}.parameter = data{idat}.parameter + 2*pattern(:,2);
  end
end

% do stats - montecarlo
cfg = [];
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05; 
cfg.numrandomization = 'all';
cfg.design = [ones(1,5) ones(1,5).*2; 1:5 1:5;];
cfg.ivar   = 1;
cfg.uvar   = 2;

% - data
cfg.parameter = 'parameter';
cfg.channel = data{1}.label;
cfg.connectivity = [...
  false true  false false false;...
  true  false false false false;...
  false false false true  false;...
  false false true  false false;...
  false false false false false;...
  ]; % neighbours: ch1 - ch2, ch3 - ch4

cfg.dim = size(data{1}.(cfg.parameter));
cfg.dimord = data{1}.dimord;
dat = cell2mat(cellfun(@(x) x.(cfg.parameter)(:), data, 'UniformOutput', false));

% - run
cfg.correctm = 'cluster';
stat = ft_statistics_montecarlo(cfg,dat,cfg.design);
cfg.correctm = 'tfce';
stat = ft_statistics_montecarlo(cfg,dat,cfg.design);


%%%%%%%%%%%%%5
% some other code


pattern = [...
    1 1 1 1 0 1 1 1 1 1;...
    1 0 0 0 0 0 0 1 0 0;...
    1 1 1 0 0 0 0 1 0 0;...
    1 0 0 0 0 0 0 1 0 0;...
    1 0 0 0 0 0 0 1 0 0;...
    ]; % FT

%% 1D (per channel)
% make fake dataset
data = cell(1,10);
for idat = 1:10
  data{idat}.label = {'ch1' 'ch2' 'ch3' 'ch4' 'ch5'}';
  data{idat}.dimord = 'chan_freq';
  data{idat}.parameter = rand(5,10);
  if idat <= 5
    data{idat}.parameter = data{idat}.parameter + 2*pattern;
  end
end
stat = run_stat(data);

% you should see the pattern in
stat.prob < 0.05

%% 2D (per channel)
% make fake dataset
data = cell(1,10);
for idat = 1:10
  data{idat}.label = {'ch1' 'ch2' 'ch3' 'ch4' 'ch5'}';
  data{idat}.dimord = 'chan_freq_time';
  data{idat}.parameter = rand(5,5,10);
  if idat <= 5
    for ch = [1 3 5]
      data{idat}.parameter(ch,:,:) = data{idat}.parameter(ch,:,:) + 2*shiftdim(pattern,-1);
    end
  end
end
stat = run_stat(data);

% you should see the pattern in
squeeze(mean(stat.prob([1 3 5],:,:),1)) < 0.05

%% 3D (per channel)
% make fake dataset
data = cell(1,10);
for idat = 1:10
  data{idat}.label = {'ch1' 'ch2' 'ch3' 'ch4' 'ch5'}';
  data{idat}.dimord = 'chan_freqlow_freqhigh_phase';
  data{idat}.parameter = rand(5,5,10,10);
  if idat <= 5
    for ch = [1 3 5]
      for ph = [2 5 8]
        data{idat}.parameter(ch,:,:,ph) = data{idat}.parameter(ch,:,:,ph) + 2*shiftdim(pattern,-1);
      end
    end
  end
end
stat = run_stat(data);

% you should see the pattern in
squeeze(mean(stat.prob([1 3 5],:,:,[2 5 8]),[1 4])) < 0.05

function stat = run_stat(data)
% do stats - montecarlo
cfg = [];
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05; 
cfg.correctm    = 'cluster'; 
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'parametric';
cfg.numrandomization = 500;
cfg.design = [ones(1,5) ones(1,5).*2; 1:5 1:5;];
cfg.ivar   = 1;
cfg.uvar   = 2;

% - data
cfg.parameter = 'parameter';
cfg.channel = data{1}.label;
cfg.connectivity = [...
  false true  false false false;...
  true  false false false false;...
  false false false true  false;...
  false false true  false false;...
  false false false false false;...
  ]; % neighbours: ch1 - ch2, ch3 - ch4

cfg.dim = size(data{1}.(cfg.parameter));
cfg.dimord = data{1}.dimord;
dat = cell2mat(cellfun(@(x) x.(cfg.parameter)(:), data, 'UniformOutput', false));

% - run
stat = ft_statistics_montecarlo(cfg,dat,cfg.design);

for fn = fieldnames(stat)'
    if size(stat.(fn{1}),1) == prod(cfg.dim)
        stat.(fn{1}) = reshape(stat.(fn{1}), cfg.dim);
    end
end
