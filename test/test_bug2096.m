function test_bug2096

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug2096
% TEST ft_sourcewrite

load(dccnpath('/home/common/matlab/fieldtrip/data/test/CP10168_4DEXP_3-Restin_BNN_V1_MEG_icaimagcoh_freq3.mat'));

cfg = [];
cfg.filename = [tempname '.cii'];
ft_sourcewrite(cfg, source);

