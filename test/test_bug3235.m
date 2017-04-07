function test_bug3235

% MEM 2gb
% WALLTIME 00:20:00

% TEST ft_volumereslice ft_sourceplot volumepermute volumeflip

%%

mri = [];
mri.dim = [11 21 31];
mri.unit = 'cm';
mri.cfg.whatever = 'something'; % for provenance tracking
mri.anatomy    = zeros(mri.dim);
mri.functional = zeros(mri.dim);

% this is a nicely aligned volume
mri.transform = [
  1 0 0 -5
  0 1 0 -10
  0 0 1 -15
  0 0 0 1
  ];

% flip it in all three directions
flip = [
  -1  0  0  0
  0 -1  0  0
  0  0 -1  0
  0  0  0  1
  ];
mri.transform = flip * mri.transform;

[X, Y, Z] = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));
pos = ft_warp_apply(mri.transform, [X(:) Y(:) Z(:)]);
mri.functional(:) = sqrt(sum(pos.^2,2)) .* prod(sign(pos),2);

cfg = [];
% cfg.anaparameter = 'anatomy';
cfg.funparameter = 'functional';
cfg.locationcoordinates = 'voxel';
cfg.location = mri.dim;
ft_sourceplot(cfg, mri);

%%

cfg = [];
cfg.method = [];
resliced = ft_volumereslice(cfg, mri);

assert(isfield(resliced, 'anatomy'));
assert(isfield(resliced, 'functional'));

% check the provenance
assert(isequal(resliced.cfg.previous, mri.cfg));
assert(~isempty(resliced.cfg.method));

cfg = [];
% cfg.anaparameter = 'anatomy';
cfg.funparameter = 'functional';
cfg.locationcoordinates = 'voxel';
cfg.location = resliced.dim;
ft_sourceplot(cfg, resliced);

%%

cfg = [];
cfg.method = 'flip';
resliced = ft_volumereslice(cfg, mri);

pos1 = ft_warp_apply(resliced.transform, [1 1 1]);      % voxel [1 1 1];
pos2 = ft_warp_apply(resliced.transform, resliced.dim); % voxel [nx ny nz];
assert(all(pos2>pos1));

cfg = [];
% cfg.anaparameter = 'anatomy';
cfg.funparameter = 'functional';
cfg.locationcoordinates = 'voxel';
cfg.location = mri.dim;
ft_sourceplot(cfg, resliced);

%%
% Arjen reported that with subject IR32 it fails
% it is a 512?666?85 CT scan, which is ~30 degrees tilted forward and which has very unisotropic voxels

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3235/IR32_CT_coreg.nii');

ct = ft_read_mri(filename);

cfg = [];
cfg.anaparameter = 'anatomy';
ft_sourceplot(cfg, ct);

cfg = [];
cfg.method = 'flip';
resliced = ft_volumereslice(cfg, ct);

cfg = [];
cfg.anaparameter = 'anatomy';
ft_sourceplot(cfg, resliced);

%% the provenance handling changed a bit, hence I also want to check that

cfg = [];
cfg.comment = filename;
ct = ft_annotate(cfg, ct);

cfg = [];
cfg.method = 'flip';
cfg.downsample = 2;
resliced = ft_volumereslice(cfg, ct);

cfg = [];
ft_analysispipeline(cfg, resliced)

% ft_volumedownsample (called inside ft_volumereslice) should not be visible in the history
assert(strcmp(resliced.cfg.previous.version.name, which('ft_annotate')))
assert(strcmp(resliced.cfg.version.name, which('ft_volumereslice')))

%%
% check the old implicit default

cfg = [];
resliced = ft_volumereslice(cfg, mri);

