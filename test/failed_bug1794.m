function failed_bug1794

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1794
% TEST ft_prepare_sourcemodel ft_convert_units ft_prepare_leadfield

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/vol.mat'))
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/dataFIC'))
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri'))

cfg = [];
cfg.grid.resolution = 1;
cfg.mri = segmentedmri;
grid = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.grid = grid;
cfg.vol = vol;
cfg.channel = dataFIC.label;
cfg.grad = dataFIC.grad;

cfgmm = cfg;
cfgmm.grid = ft_convert_units(cfg.grid,'mm');
cfgmm.vol = ft_convert_units(cfg.vol,'mm');
cfgmm.grad = ft_convert_units(cfg.grad,'mm');

cfgcm = cfg;
cfgcm.grid = ft_convert_units(cfg.grid,'cm');
cfgcm.vol = ft_convert_units(cfg.vol,'cm');
cfgcm.grad = ft_convert_units(cfg.grad,'cm');

cfgmix1 = cfg;
cfgmix1.grid = ft_convert_units(cfg.grid,'cm');
cfgmix1.vol = ft_convert_units(cfg.vol,'mm');
cfgmix1.grad = ft_convert_units(cfg.grad,'mm');

cfgmix2 = cfg;
cfgmix2.grid = ft_convert_units(cfg.grid,'mm');
cfgmix2.vol = ft_convert_units(cfg.vol,'cm');
cfgmix2.grad = ft_convert_units(cfg.grad,'cm');

cfgmix3 = cfg;
cfgmix3.grid = ft_convert_units(cfg.grid,'mm');
cfgmix3.vol = ft_convert_units(cfg.vol,'mm');
cfgmix3.grad = ft_convert_units(cfg.grad,'cm');

% % first using current (today 20130123)
% cd(dccnpath('/home/common/matlab/fieldtrip'))
% restoredefaultpath;
% ft_defaults;
gridLFmmcur = ft_prepare_leadfield(cfgmm);
gridLFcmcur = ft_prepare_leadfield(cfgcm);
gridLFmmcur.leadfield{find(gridLFmmcur.inside, 1)}(1,:)
gridLFcmcur.leadfield{find(gridLFcmcur.inside, 1)}(1,:)

% save grid1794.mat gridLF*
gridLFmix1 = ft_prepare_leadfield(cfgmix1);
gridLFmix2 = ft_prepare_leadfield(cfgmix2);
gridLFmix3 = ft_prepare_leadfield(cfgmix3);

gridLFcmcur.leadfield{find(gridLFcmcur.inside, 1)}(1,:)
gridLFmmcur.leadfield{find(gridLFmmcur.inside, 1)}(1,:)
gridLFmix1.leadfield{find(gridLFcmcur.inside, 1)}(1,:)
gridLFmix2.leadfield{find(gridLFmmcur.inside, 1)}(1,:)
gridLFmix3.leadfield{find(gridLFmmcur.inside, 1)}(1,:)

assert(isequaln(rmfield(gridLFmmcur,'cfg'),rmfield(gridLFcmcur,'cfg')))
assert(isequaln(rmfield(gridLFmix1,'cfg'),rmfield(gridLFcmcur,'cfg')))
assert(isequaln(rmfield(gridLFmix2,'cfg'),rmfield(gridLFcmcur,'cfg')))
assert(isequaln(rmfield(gridLFmix3,'cfg'),rmfield(gridLFcmcur,'cfg')))

return
%% User-provided cfg data
% this is not run in the dashboard, since file is large (3MB)

load bug1794cfg
cfgmm = cfg;
cfgmm.grid = ft_convert_units(cfg.grid,'mm');
cfgmm.vol = ft_convert_units(cfg.vol,'mm');
cfgmm.grad = ft_convert_units(cfg.grad,'mm');

cfgcm = cfg;
cfgcm.grid = ft_convert_units(cfg.grid,'cm');
cfgcm.vol = ft_convert_units(cfg.vol,'cm');
cfgcm.grad = ft_convert_units(cfg.grad,'cm');

cfgmix1 = cfg;
cfgmix1.grid = ft_convert_units(cfg.grid,'cm');
cfgmix1.vol = ft_convert_units(cfg.vol,'mm');
cfgmix1.grad = ft_convert_units(cfg.grad,'mm');

cfgmix2 = cfg;
cfgmix2.grid = ft_convert_units(cfg.grid,'mm');
cfgmix2.vol = ft_convert_units(cfg.vol,'cm');
cfgmix2.grad = ft_convert_units(cfg.grad,'cm');

cfgmix3 = cfg;
cfgmix3.grid = ft_convert_units(cfg.grid,'mm');
cfgmix3.vol = ft_convert_units(cfg.vol,'mm');
cfgmix3.grad = ft_convert_units(cfg.grad,'cm');

% % first using current (today 20130123)
% cd(dccnpath('/home/common/matlab/fieldtrip'))
% restoredefaultpath;
% ft_defaults;
gridLFmmcur = ft_prepare_leadfield(cfgmm);
gridLFcmcur = ft_prepare_leadfield(cfgcm);
gridLFmmcur.leadfield{find(gridLFmmcur.inside, 1)}(1,:)
gridLFcmcur.leadfield{find(gridLFcmcur.inside, 1)}(1,:)

% save  grid1794.mat gridLF*
gridLFmix1 = ft_prepare_leadfield(cfgmix1);
gridLFmix2 = ft_prepare_leadfield(cfgmix2);
gridLFmix3 = ft_prepare_leadfield(cfgmix3);

assert(isequaln(rmfield(gridLFmmcur,'cfg'),rmfield(gridLFcmcur,'cfg')))
assert(isequaln(rmfield(gridLFmix1,'cfg'),rmfield(gridLFcmcur,'cfg')))
assert(isequaln(rmfield(gridLFmix2,'cfg'),rmfield(gridLFcmcur,'cfg')))
assert(isequaln(rmfield(gridLFmix3,'cfg'),rmfield(gridLFcmcur,'cfg')))



%% Testing different versions

% cd(dccnpath('/home/common/matlab/fieldtrip-20120630'))
% restoredefaultpath;
% ft_defaults;
% gridLFmmjune = ft_prepare_leadfield(cfgmm);
% gridLFcmjune = ft_prepare_leadfield(cfgcm);
% 
% 
% % Here Add fieltrip-20120426 directory
% 
% cd(dccnpath('/home/common/matlab/fieldtrip'))
% restoredefaultpath;
% clear ft_defaults
% ft_defaults;
% gridLFcmcur = ft_prepare_leadfield(cfg);
% 
% 
% gridLFcmcur.leadfield{find(gridLFcmcur.inside, 1)}(1,:)
% gridLFmmcur.leadfield{find(gridLFmmcur.inside, 1)}(1,:)
% gridLFcmjune.leadfield{gridLFcmjune.inside(1)}(1,:)
% gridLFmmjune.leadfield{gridLFmmjune.inside(1)}(1,:)
% 
% 
% 
% % % Here Add fieltrip-20120426 directory
% % ft_defaults;
% % gridLFold = ft_prepare_leadfield(cfg);
% % restoredefaultpath;
% % % Here Add fieltrip-20121025 directory
% % ft_defaults;
% % cfg.grid = ft_convert_units(cfg.grid,'cm');
% % cfg.vol = ft_convert_units(cfg.vol,'cm');
% % cfg.grad = ft_convert_units(cfg.grad,'cm');
% % gridLFnew =  ft_prepare_leadfield(cfg);
% 
% 
% 
