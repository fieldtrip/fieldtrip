function test_bug2096

% WALLTIME 00:03:01

% TEST test_bug2096
% TEST ft_sourcewrite

load(dccnfilename('/home/common/matlab/fieldtrip/data/test/CP10168_4DEXP_3-Restin_BNN_V1_MEG_icaimagcoh_freq3.mat'));

cfg = [];
cfg.filename = [tempname '.cii'];
ft_sourcewrite(cfg, source);

