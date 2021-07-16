function test_ft_appendsource

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_appendsource

gridsize = 1320;
nsamples = 20;
source = [];
source.pos = randn(gridsize,3);
source.time = 1:nsamples;
source.pow = randn(gridsize,nsamples);
source.inside = 1:gridsize/2;

% test append on rpt (only case that works yet)
source2 = source;
source2.pow = randn(gridsize,nsamples);
cfg = [];
cfg.parameter = 'pow';
cfg.appenddim = 'rpt';
sourceout = ft_appendsource(cfg, source, source2);

return; %following cases do not work yet (cf. issue #1833)
%%
% test case of different positions, same time frames
source2 = source;
source2.pow = randn(gridsize,nsamples);
source2.pos = randn(gridsize,3);
cfg = [];
cfg.parameter = 'pow';
sourceout = ft_appendsource(cfg, source, source2);

% test case of different time frames, same positions
source2 = source;
source2.pow = randn(gridsize,nsamples);
source2.time = nsamples+1:nsamples*2;
cfg = [];
cfg.parameter = 'pow';
cfg.appenddim = 'rpt';
sourceout = ft_appendsource(cfg, source, source2);
