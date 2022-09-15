function inspect_bug3375

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_realtime_headlocalizer

%%

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3375'));

mri = ft_read_mri('ctf/T1/o20150923_103128t1mpragesagiso1mmwselnfps004a1001.img');

cfg = [];
cfg.method = 'interactive';
cfg.viewmode = 'surface';
cfg.coordsys = 'ctf';
mri = ft_volumerealign(cfg, mri);

cfg = [];
cfg.output = 'scalp';
seg = ft_volumesegment(cfg, mri);

cfg = [];
cfg.tissue = 'scalp';
cfg.numvertices = 10000;
scalp = ft_prepare_mesh(cfg, seg);

%%

% dewar will load automatically
% polhemus points will load automatically

cfg = [];
cfg.dataset = 'ctf/muenster_A1331.ds';
cfg.bufferdata = 'first';
cfg.headshape = scalp;
ft_realtime_headlocalizer(cfg);


%%

% dewar will load automatically
% polhemus points will load automatically

cfg = [];
cfg.dataset = 'elekta/jn_multimodal_chpi_raw_sss.fif';
cfg.bufferdata = 'first';
cfg.headshape = 'elekta/strio_jn.mat';
ft_realtime_headlocalizer(cfg);
