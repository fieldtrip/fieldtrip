%% OpenMEEG for MEG from Fieldtrip demo script
%
% This script provides an example of how to compute an MEG leadfield with
% OpenMEEG in the Fieldtrip toolbox.
%
% This demo uses spherical head models and compares the OpenMEEG result
% with the solution provided by the solution from Nolte.

%% Set the radius and conductivities of each of the compartments

addpath(cd) % Make sure current folder is in the path

close all
clear
clc

% 3 Layers
r = [88 92 100];
c = [1 1/80 1];

% % 2 Layers
% r = [92 100];
% c = [1/4 1];

% % 1 Layers
% r = [100];
% c = [1];

%% Description of the spherical mesh
[pos, tri] = icosahedron42;
% [pos, tri] = icosahedron162;
% [pos, tri] = icosahedron642;

%% Create a set of magnetometers outside the outer surface
sens.pos = max(r) * pos * 1.2;
sens.ori = pos;
sens.label = {};
nsens = size(sens.pos, 1);
for ii=1:nsens
    sens.label{ii} = sprintf('vertex%03d', ii);
end

%% Set the position of the probe dipole
pos = [0 0 70];

%% Create a BEM volume conduction model
vol = [];
vol1 = [];
for ii=1:length(r)
    vol.bnd(ii).pos = pos * r(ii);
    vol.bnd(ii).tri = tri;
    if (ii==1);
        vol1.bnd(ii).pos = pos * r(ii);
        vol1.bnd(ii).tri = tri;
    end
end
vol.cond = c;
vol1.cond = c(1);

%% choose MEG implementation (Nolte, OpenMEEG)

% Compute the BEM
cfg.method = 'openmeeg';
vol = ft_prepare_bemmodel(cfg, vol);

cfg.vol = vol;
cfg.grid.pos = pos;
cfg.grad = sens;
cfg.reducerank = 'no';
grid = ft_prepare_leadfield(cfg);
lf_openmeeg = grid.leadfield{1};

% choose MEG Nolte
clear cfg;
cfg.method = 'singleshell';
cfg.grid.pos = pos;
cfg.grad = sens;
vol1.type = 'singleshell';
[vol1,sens] = ft_prepare_vol_sens(vol1, sens);
cfg.vol = vol1;
cfg.reducerank = 'no';
grid = ft_prepare_leadfield(cfg);
lf_singleshell = grid.leadfield{1};

%% Plot both OpenMEEG and analytic leadfield for visual inspection
figure
hold on
plot(lf_openmeeg(:,1),'bx-','linewidth',2)
plot(lf_singleshell(:,1),'r--','linewidth',2)
hold off
legend({'OpenMEEG' 'Nolte'})
