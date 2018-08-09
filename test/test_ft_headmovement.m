function test_ft_headmovement

% MEM 3000mb
% WALLTIME 00:20:00

% TEST ft_headmovement

dataset = dccnpath('/home/common/matlab/fieldtrip/data/SubjectRest.ds');

cfg = [];
cfg.dataset = dataset;
data1 = ft_headmovement(cfg);

cfg = [];
cfg.dataset = dataset;
cfg.method = 'cluster';
cfg.numclusters = 10;
data2 = cell(1,cfg.numclusters);
[data2{:}] = ft_headmovement(cfg);

cfg = [];
cfg.length = 10;
data_chunked = ft_redefinetrial(cfg, data1);

cfg = [];
cfg.dataset = dataset;
cfg.trl    = data_chunked.sampleinfo;
cfg.trl(:,3) = 0;
cfg.method = 'pertrial';
data3 = cell(1,numel(data_chunked.trial));
[data3{:}] = ft_headmovement(cfg);
