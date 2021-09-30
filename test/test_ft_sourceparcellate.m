function test_ft_sourceparcellate

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_sourceparcellate ft_read_atlas ft_sourceinterpolate

atlasfilename   = dccnpath('/home/common/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');

gridsize = 38*48*41;
nsamples = 20;
source = [];
source.pos = randn(gridsize,3);
source.time = 1:nsamples;
source.pow = randn(gridsize,nsamples);
source.inside = 1:gridsize;

atlas = ft_read_atlas(atlasfilename);
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
parcel = ft_sourceinterpolate(cfg, atlas, source);
parcel.pos = source.pos;

cfg = [];
cfg.parcellation = 'tissue';
meshout = ft_sourceparcellate(cfg, source, parcel);