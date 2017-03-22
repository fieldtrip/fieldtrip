function test_bug2185

% WALLTIME 00:20:00
% MEM 6gb

% TEST ft_sourcegrandaverage ft_selectdata ft_selectdata_new ft_datatype_source ft_math

global ft_default
ft_default = [];

%% load the end-user provided test data

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2185.mat');
load(filename);

% The content of the file is a cell-array that looks like this
%
% source_timelock_stim{1}
% ans =
%        time: [1x1500 double]
%         pos: [8196x3 double]
%      inside: [8196x1 double]
%     outside: [1x0 double]
%      method: 'average'
%         avg: [1x1 struct]
%         cfg: [1x1 struct]
% source_timelock_stim{1}.avg
% ans =
%          mom: {8196x1 cell}
%          pow: [8196x1500 double]
%     noisecov: {8196x1 cell}

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});
assert(isfield(output, 'pow'), 'missing output field');
assert(isfield(output, 'time'), 'missing output field');

cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'multiply';
cfg.scalar = 0;
output = ft_math(cfg, output);
assert(all(output.pow(:)==0), 'power should be zero');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intermezzo, do some interpolation and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source2d = output;
source2d.pow(:,:) = 0;
source2d.pow(:,1) = source2d.pos(:,3); % replace the values by something easier to visualize

% create a regular 3d-grid source structure that encompasses the cortical sheet
minpos = floor(min(source2d.pos));
maxpos = ceil(max(source2d.pos));
xgrid = minpos(1):0.5:maxpos(1);
ygrid = minpos(2):0.5:maxpos(2);
zgrid = minpos(3):0.5:maxpos(3);
dim = [length(xgrid) length(ygrid) length(zgrid)];
[X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
source3d.pos = [X(:) Y(:) Z(:)];
source3d.dim = dim;
% convert the regular 3d-grid source into a volume
volume3d = ft_checkdata(source3d, 'datatype', 'volume');

cfg = [];
cfg.parameter = 'pow';
source3d = ft_sourceinterpolate(cfg, source2d, volume3d);

cfg = [];
cfg.avgovertime = 'yes';
cfg.keeptimedim = 'no';
cfg.parameter = 'pow';
source3davg = ft_selectdata(cfg, source3d);

cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg, source3davg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of intermezzo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'yes';
output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});
assert( isfield(output, 'pow'), 'missing output field');
assert(isfield(output, 'time'), 'missing output field');

cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'multiply';
cfg.scalar = 0;
output = ft_math(cfg, output);
assert(all(output.pow(:)==0), 'power should be zero');

cfg = [];
cfg.parameter = 'mom';
cfg.keepindividual = 'no';
output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});
assert( isfield(output, 'mom'), 'missing output field');
assert( isfield(output, 'time'), 'missing output field');

cfg = [];
cfg.parameter = 'mom';
cfg.operation = 'multiply';
cfg.scalar = 0;
output = ft_math(cfg, output);
assert(all(output.mom{1}(:)==0), 'moment should be zero');

cfg = [];
cfg.parameter = 'mom';
cfg.keepindividual = 'yes';
output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});
assert( isfield(output, 'mom'), 'missing output field');
assert( isfield(output, 'time'), 'missing output field');

cfg = [];
cfg.parameter = 'mom';
cfg.operation = 'multiply';
cfg.scalar = 0;
output = ft_math(cfg, output);
assert(all(output.mom{1}(:)==0), 'moment should be zero');

cfg = [];
cfg.parameter = 'noisecov';
cfg.keepindividual = 'yes';
output = ft_sourcegrandaverage(cfg, source_timelock_stim{:});

