function inspect_bug1223

% TEST inspect_bug1223
% TEST ft_read_mri ft_sourceplot

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1223/Num1_MinDef_M_Normal_age12_num10Atlas.hdr'));

cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
