function test_ft_headmovement

% MEM 3gb
% WALLTIME 00:20:00
% DEPENDENCY ft_headmovement

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/SubjectRest.ds');

cfg = [];
cfg.dataset = dataset;
cfg.method  = 'updatesens';
data1 = ft_headmovement(cfg);

cfg = [];
cfg.dataset = dataset;
cfg.method = 'cluster';
cfg.numclusters = 10;
data2 = cell(1,cfg.numclusters);
[data2{:}] = ft_headmovement(cfg);

cfg          = [];
cfg.dataset  = dataset;
cfg.trl      = [(1:1000:100000)' (1000:1000:100000)'];
cfg.trl(:,3) = 0;
cfg.method = 'pertrial_cluster';
cfg.numclusters = 10;
data3 = cell(1,cfg.numclusters);
[data3{:}] = ft_headmovement(cfg);

cfg         = [];
cfg.dataset = dataset;
cfg.trl      = [(1:1000:100000)' (1000:1000:100000)'];
cfg.trl(:,3) = 0;
cfg.method  = 'avgoverrpt';
data4 = ft_headmovement(cfg);
