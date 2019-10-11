function test_issue1198

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_sourceanalysis ft_sourcedescriptives ft_sourceplot

%%
% http://www.fieldtriptoolbox.org/example/testing_bem_created_leadfields/

% Create a spherical volume conductor of radius 100, conductivity 1 and center [0,0,0]
vol   = [];
vol.r = 100;
vol.c = 1;
vol.o = [0 0 0];
vol.unit = 'mm';

% Create a set of electrodes on the upper half of the sphere
[X, Y, Z] = sphere(10);
pos = unique([X(:) Y(:) Z(:)], 'rows');
pos = pos(pos(:,3)>=0,:);

elec = [];
elec.elecpos = vol.r * pos;
elec.label = {};
nelec = size(pos,1);
for ii = 1:nelec
  elec.label{ii} = sprintf('vertex%03d', ii);
end

% Define positions of dipoles (along z axis)
zp = linspace(0,100,50)';
pos = [zeros(size(zp)) zeros(size(zp)) zp];

% Define the corresponding spatial grid
cfg = [];
cfg.inwardshift = 5;
cfg.sourcemodel.pos = pos;
cfg.headmodel = vol;
cfg.elec = elec;
sourcemodel = ft_prepare_sourcemodel(cfg);

% construct a dense triangulated mesh
[X, Y, Z] = sphere(100);
pos = unique([X(:) Y(:) Z(:)], 'rows');
mesh.pos = pos*100;
mesh.tri = convhulln(mesh.pos);

vertic = [50 80 100 200 2000];
for ll=1:length(vertic)
  cfg = [];
  cfg.method = 'dipoli';
  cfg.conductivity  = 1;
  cfg.numvertices = vertic(ll);
  cfg.isolatedsource = false;
  volbem{ll} = ft_prepare_headmodel(cfg, mesh);
end

% calculate BEM leadfield
for ll=1:length(vertic)
  cfg = [];
  cfg.sourcemodel = sourcemodel;
  cfg.headmodel = volbem{ll};
  cfg.elec = elec;
  leadfieldBEM{ll} = ft_prepare_leadfield(cfg);
end

% calculate theoretical single sphere leadfield
cfg = [];
cfg.grid = sourcemodel;
cfg.headmodel = vol;
cfg.elec = elec;
leadfieldSphere = ft_prepare_leadfield(cfg);


%%
% http://www.fieldtriptoolbox.org/example/compute_forward_simulated_data_and_apply_a_beamformer_scan/

% This example script shows you how to create some simulated channel-level
% MEG data with a single dipole at a specified location in the head.
% Subsequently it does a beamformer source reconstruction to localize that
% source.

% create an array with some magnetometers at 12cm distance from the origin
[X, Y, Z] = sphere(10);
pos = unique([X(:) Y(:) Z(:)], 'rows');
pos = pos(pos(:,3)>=0,:);
grad = [];
grad.coilpos = 12*pos;
grad.coilori = pos; % in the outward direction
% grad.tra = eye(length(pos)); % each coils contributes exactly to one channel
for i=1:length(pos)
  grad.label{i} = sprintf('chan%03d', i);
end

% create a spherical volume conductor with 10cm radius
vol.r = 10;
vol.o = [0 0 0];

% note that beamformer scanning will be done with a 1cm grid, so you should
% not put the dipole on a position that will not be covered by a grid
% location later
cfg = [];
cfg.headmodel = vol;
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
cfg.headmodel = vol;
cfg.grad = grad;
cfg.resolution = 1;
cfg.method = 'lcmv';
cfg.lcmv.projectnoise = 'yes'; % needed for neural activity index
source = ft_sourceanalysis(cfg, timelock);

% compute the neural activity index, i.e. projected power divided by
% projected noise
cfg = [];
cfg.powmethod = 'none'; % keep the power as estimated from the data covariance, i.e. the induced power
source_nai = ft_sourcedescriptives(cfg, source);

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'nai';
cfg.funcolorlim = [1.4 1.5];  % the voxel in the center of the volume conductor messes up the autoscaling
ft_sourceplot(cfg, source_nai);

