% function test_example_sourcerecon_meeg

% MEM 2gb
% WALLTIME 00:10:00
% DATA none

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a set of electrodes, randomly placed on the sphere
elec = [];
elec.elecpos = randn(32,3);
dum = sqrt(sum(elec.elecpos.^2,2));                  % compute the distance to the origin
elec.elecpos = elec.elecpos ./ [dum dum dum] * 0.10; % scale them to a 0.1 meter sphere
for i=1:32
  elec.label{i} = sprintf('eeg%03d', i);
end

% create a concentric 3-sphere volume conductor for the EEG, the radius is the same as for the electrodes
headmodel_eeg   = [];
headmodel_eeg.r = [0.088 0.092 0.100]; % radii of spheres
headmodel_eeg.c = [1 1/80 1];          % conductivity
headmodel_eeg.o = [0 0 0];             % center of sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a set of magnetometers, randomly placed around the sphere
grad = [];
grad.coilpos = randn(64,3);
dum = sqrt(sum(grad.coilpos.^2,2));
grad.coilpos = grad.coilpos ./ [dum dum dum] * 0.12; % scale them to a 0.12 meter sphere, shifted outward from the head surface
grad.coilori = grad.coilpos ./ [dum dum dum];        % unit length
for i=1:64
  grad.label{i} = sprintf('meg%03d', i);
end
grad.tra = eye(64,64);

% create a single-sphere volume conductor for the MEG
headmodel_meg   = [];
headmodel_meg.r = 0.1;             % radius of sphere
headmodel_meg.c = 1;                % conductivity
headmodel_meg.o = [0 0 0];          % center of sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine the EEG and MEG sensor definitions and volume conductor models
% and do a forward computation
combined_headmodel = {headmodel_eeg, headmodel_meg};
combined_sens = {elec, grad};

pos = [0 0 0.8];
mom = [1 0 0]';

[combined_headmodel{1}, combined_sens{1}] = ft_prepare_vol_sens(combined_headmodel{1}, combined_sens{1});
[combined_headmodel{2}, combined_sens{2}] = ft_prepare_vol_sens(combined_headmodel{2}, combined_sens{2});
leadfield  = ft_compute_leadfield(pos, combined_sens, combined_headmodel);

assert(size(leadfield,1)==96);
assert(size(leadfield,2)==3);

eeg = leadfield(1:32,:)*mom;
meg = leadfield(33:end,:)*mom;

figure; plot(eeg); title('eeg');
figure; plot(meg); title('meg');

ft_senstype(combined_sens)
ft_headmodeltype(combined_headmodel)

figure
ft_plot_topo3d(elec.elecpos, eeg)
ft_plot_sens(elec, 'axes', true)
ft_plot_headmodel(headmodel_eeg, 'facecolor', 'skin', 'facealpha', 0.5)

figure
ft_plot_topo3d(grad.coilpos, meg)
ft_plot_sens(grad, 'axes', true)
ft_plot_headmodel(headmodel_meg, 'facecolor', 'skin', 'facealpha', 0.5)

