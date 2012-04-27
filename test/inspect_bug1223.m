function test_bug1223

% TEST test_bug1223
% TEST ft_read_mri ft_sourceplot

cd H:/common/matlab/fieldtrip/data/test/bug1223

mri = ft_read_mri('Num1_MinDef_M_Normal_age12_num10Atlas.hdr');

cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
