function test_ft_combineplanar

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY ft_combineplanar ft_preprocessing ft_timelockanalysis ft_prepare_neighbours ft_megplanar

% template loading should be modified once template MEG data are available
% (cf.issue #1834)
subjectfilename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');

cfg            = [];
cfg.dataset    = subjectfilename;
cfg.channel    = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
data = ft_preprocessing(cfg);

cfg = [];
data = ft_timelockanalysis(cfg,data);

cfg              = [];
cfg.feedback     = 'no';
cfg.method       = 'template';
cfg.neighbours   = ft_prepare_neighbours(cfg, data);
cfg.planarmethod = 'sincos';
data = ft_megplanar(cfg, data);

cfg = [];
cfg.updatesens = 'no';
dataout = ft_combineplanar(cfg,data);

