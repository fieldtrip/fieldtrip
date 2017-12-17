function test_bug2556

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_sourceparcellate ft_checkdata

[ftver, ftpath] = ft_version;
filename = fullfile(ftpath, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii');

aal = ft_read_atlas(filename);

aal = ft_checkdata(aal, 'datatype', 'parcellation', 'parcellationstyle', 'indexed');

%% part 2, deal with pow and avg.pow confusion

source = [];
source.pos = randn(50,3);
source.avg.pow = randn(50,1);

parcellation = [];
parcellation.pos = source.pos;
parcellation.tissue = [ones(1,25)*1 ones(1,25)*2]';
parcellation.tissuelabel = {'left', 'right'};

cfg = [];
source1p = ft_sourceparcellate(cfg, source, parcellation);

source.pow = source.avg.pow;
source = rmfield(source, 'avg');

cfg = [];
cfg.parameter = 'pow';
source2p = ft_sourceparcellate(cfg, source, parcellation);

assert(isequal(source1p.pow,   source2p.pow));
assert(isequal(source1p.label, source2p.label));

%% part 3, check ft_sourcegrandaverage

cfg = [];
cfg.keepindividual = 'no';
grandavg = ft_sourcegrandaverage(cfg, source, source, source);

cfg = [];
cfg.parameter = 'pow';
source3p = ft_sourceparcellate(cfg, grandavg, parcellation);

% assert(isequal(source1p.pow,   source3p.pow));
assert(isalmostequal(source1p.pow,   source3p.pow, 'reltol', 1e-8)); % there is a tiny numerical difference
assert(isequal(source1p.label, source3p.label));





