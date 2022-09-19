function test_bug3369

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY read_neuromag_maxfilterlog

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3369'));

file0 = 'tactile_stim_raw.fif_tsss_mc.log';
file1 = 'tactile_stim_raw-1.fif_tsss_mc.log';
file2 = 'tactile_stim_raw-2.fif_tsss_mc.log';

%%

cfg = [];
cfg.dataset = file0;
cfg.viewmode = 'vertical';
ft_databrowser(cfg)

%%

cfg = [];
cfg.dataset = file0;
data0 = ft_preprocessing(cfg);
cfg.dataset = file1;
data1 = ft_preprocessing(cfg);
cfg.dataset = file2;
data2 = ft_preprocessing(cfg);

cfg = [];
cfg.keepsampleinfo = 'no'; % the sampleinfo of the three does not match
data = ft_appenddata(cfg, data0, data1, data2);

%%

cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
ft_databrowser(cfg, data);
