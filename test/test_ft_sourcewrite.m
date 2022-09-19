function test_ft_sourcewrite

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_sourcewrite

gridsize = 38*48*41;
nsamples = 20;
source = [];
source.dim = [38,48,41];
source.pos = randn(gridsize,3);
source.time = 1:nsamples;
source.pow = randn(gridsize,nsamples);
source.inside = 1:gridsize;

cfg = [];
cfg.parameter = 'pow';
cfg.filename = tempname;
cfg.filetype = 'nifti';
ft_sourcewrite(cfg, source)
delete([cfg.filename '.nii'])
