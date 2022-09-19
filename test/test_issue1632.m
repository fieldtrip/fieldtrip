function test_issue1632

% MEM 4gb
% WALLTIME 00:10:00
% DEPENDENCY ft_compute_leadfield

% https://github.com/fieldtrip/fieldtrip/issues/1632
% This script will demonstrate that with ft_prepare leadfields, the
% leadfields are calculated for every coil and not every sensor.
% Subsequently, it will show that ft_sourceanalysis lines 437-443 truncate
% the leadfield matrix, removing all lines beyond the number of sensors.

% The mat-file includes:
% data (emptied except for grad and label) from ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/beamformer_lcmv/
% headmodel created with ft_prepare_mesh and ft_prepare_headmodel with
% method 'openmeeg', based on segmentedmri.mat from ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/beamformer_lcmv/
% Author: B Knipscheer 13-01-2021

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1632.mat'));

cfg                  = [];
cfg.grad             = data.grad;
cfg.headmodel        = headmodel;
cfg.sourcemodel.pos  = [0 0 0.07];
cfg.channel          = 'MEG';
cfg.normalize        = 'column';
sourcemodel          = ft_prepare_leadfield(cfg);

fprintf('Number of good MEG sensors = %d\n', size(data.label,1));  % There are 274 MEG sensors
fprintf('Number of coils matching the MEG sensors = %d\n', size(data.grad.tra,2)); % There are 548 MEG coils for the 274 MEG sensors.
fprintf('Number of leadfield rows for point 1 before ft_sourceanalysis = %d\n', size(sourcemodel.leadfield{1},1));

assert(isequal(size(sourcemodel.leadfield{1},1), numel(data.label)));
