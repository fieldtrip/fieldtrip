function test_bug1870

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_ft_megplanar ft_megplanar ft_datatype_sens ft_compute_leadfield

% this test is basically a small section of test_ft_megplanar
% the input data is consistent, but along the way the grad structure gets screwed up

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1870.mat'));
dataP = ft_megplanar(cfg, data);

% although it was not reported, also check on the combining of the planar channels
cfg = [];
dataC = ft_combineplanar(cfg, dataP);

% do some sanity checks
assert(isfield(dataP.grad, 'tra'));
ft_datatype_sens(dataP.grad);
ft_datatype_sens(dataC.grad);

