function test_bug2316

% MEM 1gb
% TEST test_bug2316
% TEST ft_databrowser

% this function currently serves as a placeholder to reproduce bug2316,
% i.e. it crashes Matlab (2011a, but possibly other versions as well) on
% a PC running Windows 7. on MacOS it seems fine. 

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2316.mat');
load(filename);

cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.viewmode = 'component';
ft_databrowser(cfg, datacomp);
