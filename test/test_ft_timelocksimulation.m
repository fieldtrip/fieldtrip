function test_ft_timelocksimulation

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_timelocksimulation
% DATA no

%%
cfg = [];
dataout = ft_timelocksimulation(cfg);

%%
% the following tests the functionality of setting the random seed (or not)
cfg = [];
cfg.s1.numcycli = 1;
cfg.s1.ampl     = 1.0;
cfg.s2.numcycli = 2;
cfg.s2.ampl     = 0.7;
cfg.s3.numcycli = 4;
cfg.s3.ampl     = 0.2;
cfg.noise.ampl  = 0.1;
data1 = ft_timelocksimulation(cfg);

cfg = [];
cfg.s1.numcycli = 1;
cfg.s1.ampl     = 1.0;
cfg.s2.numcycli = 2;
cfg.s2.ampl     = 0.7;
cfg.s3.numcycli = 4;
cfg.s3.ampl     = 0.2;
cfg.noise.ampl  = 0.1;
data2 = ft_timelocksimulation(cfg);

% they should have different noise
assert(~isequal(data1.trial, data2.trial));

% now use the same random seed twice, should yield the same results
cfg = [];
cfg.s1.numcycli = 1;
cfg.s1.ampl     = 1.0;
cfg.s2.numcycli = 2;
cfg.s2.ampl     = 0.7;
cfg.s3.numcycli = 4;
cfg.s3.ampl     = 0.2;
cfg.noise.ampl  = 0.1;
cfg.randomseed  = 123;
data1 = ft_timelocksimulation(cfg);

cfg = [];
cfg.s1.numcycli = 1;
cfg.s1.ampl     = 1.0;
cfg.s2.numcycli = 2;
cfg.s2.ampl     = 0.7;
cfg.s3.numcycli = 4;
cfg.s3.ampl     = 0.2;
cfg.noise.ampl  = 0.1;
cfg.randomseed  = 123;
data2 = ft_timelocksimulation(cfg);

% they should have the same noise
assert(isequal(data1.trial, data2.trial));