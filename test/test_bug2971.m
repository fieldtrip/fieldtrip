function test_bug2971

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_sourcestatistics ft_selectdata

% should be the same for all source structures
pos = randn(10,3);

source = {};
for i=1:10
  source{i} = [];
  source{i}.dim  = [1 2 5];
  source{i}.pos  = pos;
  source{i}.time = 1:6;
  source{i}.freq = 1:10;
  source{i}.stat = randn(10,1);
  source{i}.statdimord = 'pos';
  source{i}.avg.pow = randn(10,1);
end

% ft_selectdata/getdimord has difficulties with source.stat, which is either 'freq' or 'pos'
% it does not apply to source.avg.pow, because getdimord knows the most likely dimord of the pow field

cfg = [];
cfg.design = [1 1 1 1 1 2 2 2 2 2];
cfg.ivar = 1;
cfg.parameter = 'stat';
cfg.method = 'analytic';
cfg.statistic = 'indepsamplesT';
stat = ft_sourcestatistics(cfg, source{:});

