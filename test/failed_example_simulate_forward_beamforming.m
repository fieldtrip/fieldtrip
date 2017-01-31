function failed_example_simulate_forward_beamforming

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_example_simulate_forward_beamforming
% TEST ft_dipolesimulation ft_timelockanalysis ft_sourceanalysis

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% This example script shows you how to create some simulated channel-level
% MEG data with a single dipole at a specified location in the head.
% Subsequently it does a beamformer source reconstruction to localize that
% source.

% create a gradiometer array with magnetometers at 12cm distance from the origin
[pnt, tri] = icosahedron162;
pnt = pnt(pnt(:,3)>=0,:);
grad.pnt = 12*pnt;
grad.ori = pnt;
for i=1:length(pnt)
  grad.label{i} = sprintf('chan%03d', i);
end

% create a spherical volume conductor with 10cm radius
vol.r = 10;
vol.o = [0 0 0];

% note that beamformer scanning will be done with a 1cm grid, so you should
% not put the dipole on a position that will not be covered by a grid
% location later
cfg = [];
cfg.vol = vol;
cfg.grad = grad;
cfg.dip.pos = [0 0 4];    % you can vary the location, here the dipole is along the z-axis
cfg.dip.mom = [1 0 0]';   % the dipole points along the x-axis
cfg.relnoise = 10;
cfg.ntrials = 20;
data = ft_dipolesimulation(cfg);

% compute the data covariance matrix, which will capture the activity of
% the simulated dipole
cfg = [];
cfg.covariance = 'yes';
timelock = ft_timelockanalysis(cfg, data);

% do the beamformer source reconstuction on a 1 cm grid
cfg = [];
cfg.vol = vol;
cfg.grad = grad;
cfg.grad.unit='cm'; %error otherwise
cfg.resolution = 1;
cfg.method = 'lcmv';
cfg.projectnoise='yes'; %needed for neural activity index
source = ft_sourceanalysis(cfg, timelock);

% compute the neural activity index, i.e. projected power divided by
% projected noise
cfg = [];
cfg.powmethod = 'none'; % keep the power as estimated from the data covariance, i.e. the induced power
source = ft_sourcedescriptives(cfg, source);
source.avg.nai=source.avg.pow./source.avg.noise;  %neural activity index calculation

source = rmfield(source,'time');

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'avg.nai';
cfg.funcolorlim = [1.5 2];  % the voxel in the center of the volume conductor messes up the autoscaling
ft_sourceplot(cfg, source);
