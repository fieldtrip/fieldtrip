%% OpenMEEG for EEG from Fieldtrip demo script
% This script provides an example of how to compute an EEG leadfield with OpenMEEG in the Fieldtrip toolbox.
%
% This demo uses spherical head models and compares the OpenMEEG result with the analytical solution.

%% Set the radius and conductivities of each of the compartments

addpath(cd) % Make sure current folder is in the path

close all
clear
clc

% % 4 Layers
% r = [85 88 92 100];
% c = [1 1/20 1/80 1];

% 3 Layers
r = [88 92 100];
c = [1 1/80 1];

% % 2 Layers
% r = [100 92];
% c = [1 1/4];
% 
% % 1 Layers
% r = [100];
% c = [1];

%% Description of the spherical mesh
[pos, tri] = mesh_sphere(42);
% [pos, tri] = mesh_sphere(162);
% [pos, tri] = mesh_sphere(642);

%% Create a set of electrodes on the outer surface
sens.elecpos = max(r) * pos;
sens.label = {};
nsens = size(sens.elecpos,1);
for ii=1:nsens
    sens.label{ii} = sprintf('vertex%03d', ii);
end

%% Set the position of the probe dipole
dip_pos = [0 0 70];

%% Create a BEM volume conduction model
vol = [];
for ii=1:length(r)
    vol.bnd(ii).pos = pos * r(ii);
    vol.bnd(ii).tri = fliplr(tri); % pointing inwards!!!
end

%% Compute the BEM

% choose BEM implementation (OpenMEEG, bemcp or dipoli)
% cfg=[];
% cfg.method = 'openmeeg';
% vol = ft_prepare_bemmodel(cfg, vol);

cfg=[];
cfg.method = 'openmeeg';
cfg.conductivity = c;
vol = ft_prepare_headmodel(cfg, vol);

cfg.headmodel = vol;
cfg.grid.pos = dip_pos;
cfg.elec = sens;
grid = ft_prepare_leadfield(cfg);

lf_openmeeg = grid.leadfield{1};

%% Plot result
bnd = struct('pos', pos, 'tri', tri);
figure; ft_plot_mesh(bnd, 'vertexcolor', lf_openmeeg(:,1))
figure; ft_plot_mesh(bnd, 'vertexcolor', lf_openmeeg(:,2))
figure; ft_plot_mesh(bnd, 'vertexcolor', lf_openmeeg(:,3))

%% Compute the analytic leadfield

vol_sphere.r = r;
vol_sphere.cond = c;

lf_sphere = ft_compute_leadfield(dip_pos, sens, vol_sphere);

%% Evaluate the quality of the result using RDM and MAG
rdms = zeros(1,size(lf_openmeeg,2));
for ii=1:size(lf_openmeeg,2)
    rdms(ii) = norm(lf_openmeeg(:,ii)/norm(lf_openmeeg(:,ii)) - lf_sphere(:,ii)/norm(lf_sphere(:,ii)));
end
mags = sqrt(sum(lf_openmeeg.^2))./sqrt(sum(lf_sphere.^2));
disp(['RDMs: ',num2str(rdms)]);
disp(['MAGs: ',num2str(mags)]);

%% Plot both OpenMEEG and analytic leadfield for visual inspection
figure
plot([lf_openmeeg(:,1), lf_sphere(:,1)])
legend({'OpenMEEG' 'Analytic'})
