function test_bug2539

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_datatype ft_checkdata

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2539.mat'))

cfg=[];
cfg.component          = 1:16;
cfg.layout             = layout;
cfg.zlim               = 'maxabs';
cfg.marker             = 'on';

ft_topoplotIC(cfg, weights);

assert(ft_datatype(weights, 'comp'))
assert(~ft_datatype(weights, 'raw'))
assert(~ft_datatype(weights, 'timelock'))
assert(~ft_datatype(weights, 'freq'))
