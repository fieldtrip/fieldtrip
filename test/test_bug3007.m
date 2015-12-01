function test_bug3007

% WALLTIME 00:10:00
% MEM 1gb

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3007'))

load cfg
load eye_movs

[cfgout, movement] = ft_detect_movement(cfg, eye_movs);

assert(isstruct(cfgout));
assert(~isstruct(movement)); % should be Nx3 array
