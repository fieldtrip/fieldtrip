function test_ft_sourceparcellate

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_sourceparcellate ft_read_atlas ft_sourceinterpolate

% see also https://github.com/fieldtrip/fieldtrip/issues/1753

%%

atlasfilename   = dccnpath('/home/common/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');

npos = 38*48*41;
ntime = 20;

source = [];
source.pos = 200*rand(npos,3) - 100;
source.unit = 'mm';
source.time = 1:ntime;
source.pow = randn(npos,ntime);
source.inside = 1:npos;

atlas = ft_read_atlas(atlasfilename);

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
parcellation = ft_sourceinterpolate(cfg, atlas, source);

cfg = [];
cfg.parcellation = 'tissue';
result = ft_sourceparcellate(cfg, source, parcellation);

%%

x = -5:5;
y = -5:5;
z = -5:5;

[X, Y, Z] = ndgrid(x, y, z);
dim = [length(x) length(y) length(z)];

source = [];
source.dim = dim;
source.pos = [X(:) Y(:) Z(:)];
source.time = linspace(-0.1, 0.5, 15);
source.pow   = randn(prod(dim), 1);
for i=1:prod(dim)
  source.eta{i} = source.pow(i);
end
source.timecourse = randn(prod(dim), 15);
source.moment     = cell(prod(dim), 1);
for i=1:prod(dim)
  source.moment{i} = source.timecourse(i,:);
end
source.powdimord = 'pos';
source.etadimord = '{pos}';
source.timecoursedimord = 'pos_time';
source.momentdimord = '{pos}_ori_time';

atlas = [];
atlas.dim = dim;
atlas.pos = [X(:) Y(:) Z(:)];
atlas.tissue = zeros(prod(dim),1);
roi1 = atlas.pos(:,1)<0;
roi2 = atlas.pos(:,1)>0;
atlas.tissue(roi1) = 1;
atlas.tissue(roi2) = 2;
atlas.tissuelabel{1} = 'left';
atlas.tissuelabel{2} = 'right';

%%

cfg = [];
cfg.method = 'mean';

cfg.parameter = 'pow';
parcel1 = ft_sourceparcellate(cfg, source, atlas);
manual = mean(source.(cfg.parameter)(roi1, :), 1);
assert(isequal(parcel1.pow(1,:), manual));

cfg.parameter = 'eta';
parcel2 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel1.pow(:), parcel2.eta(:), 'abstol', 100*eps));

cfg.parameter = 'timecourse';
parcel3 = ft_sourceparcellate(cfg, source, atlas);
manual = mean(source.(cfg.parameter)(roi1, :), 1);
assert(isequal(parcel3.timecourse(1,:), manual));

cfg.parameter = 'moment';
parcel4 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel3.timecourse(:), parcel4.moment(:), 'abstol', 100*eps));

%%

cfg = [];
cfg.method = 'median';

cfg.parameter = 'pow';
parcel1 = ft_sourceparcellate(cfg, source, atlas);
manual = median(source.(cfg.parameter)(roi1, :), 1);
assert(isequal(parcel1.pow(1,:), manual));

cfg.parameter = 'eta';
parcel2 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel1.pow(:), parcel2.eta(:), 'abstol', 100*eps));

cfg.parameter = 'timecourse';
parcel3 = ft_sourceparcellate(cfg, source, atlas);
manual = median(source.(cfg.parameter)(roi1, :), 1);
assert(isequal(parcel3.timecourse(1,:), manual));

cfg.parameter = 'moment';
parcel4 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel3.timecourse(:), parcel4.moment(:), 'abstol', 100*eps));

%%

cfg = [];
cfg.method = 'min';

cfg.parameter = 'pow';
parcel1 = ft_sourceparcellate(cfg, source, atlas);
manual = min(source.(cfg.parameter)(roi1, :), [], 1);
assert(isequal(parcel1.pow(1,:), manual));

cfg.parameter = 'eta';
parcel2 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel1.pow(:), parcel2.eta(:), 'abstol', 100*eps));

cfg.parameter = 'timecourse';
parcel3 = ft_sourceparcellate(cfg, source, atlas);
manual = min(source.(cfg.parameter)(roi1, :), [], 1);
assert(isequal(parcel3.timecourse(1,:), manual));

cfg.parameter = 'moment';
parcel4 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel3.timecourse(:), parcel4.moment(:), 'abstol', 100*eps));

%%

cfg = [];
cfg.method = 'max';

cfg.parameter = 'pow';
parcel1 = ft_sourceparcellate(cfg, source, atlas);
manual = max(source.(cfg.parameter)(roi1, :), [], 1);
assert(isequal(parcel1.pow(1,:), manual));

cfg.parameter = 'eta';
parcel2 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel1.pow(:), parcel2.eta(:), 'abstol', 100*eps));

cfg.parameter = 'timecourse';
parcel3 = ft_sourceparcellate(cfg, source, atlas);
manual = max(source.(cfg.parameter)(roi1, :), [], 1);
assert(isequal(parcel3.timecourse(1,:), manual));

cfg.parameter = 'moment';
parcel4 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel3.timecourse(:), parcel4.moment(:), 'abstol', 100*eps));

%%

cfg = [];
cfg.method = 'std';

cfg.parameter = 'pow';
parcel1 = ft_sourceparcellate(cfg, source, atlas);
manual = std(source.(cfg.parameter)(roi1, :), 0, 1);
assert(isequal(parcel1.pow(1,:), manual));

cfg.parameter = 'eta';
parcel2 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel1.pow(:), parcel2.eta(:), 'abstol', 100*eps));

cfg.parameter = 'timecourse';
parcel3 = ft_sourceparcellate(cfg, source, atlas);
manual = std(source.(cfg.parameter)(roi1, :), 0, 1);
assert(isequal(parcel3.timecourse(1,:), manual));

cfg.parameter = 'moment';
parcel4 = ft_sourceparcellate(cfg, source, atlas);
assert(isalmostequal(parcel3.timecourse(:), parcel4.moment(:), 'abstol', 100*eps));
