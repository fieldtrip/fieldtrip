function test_example_custom_forwinv

% MEM 4gb
% WALLTIME 00:20:00

%
%% Use your own forward leadfield model in an inverse beamformer computation
%
% You can use externally precomputed leadfield matrices for each dipole location on which you want to compute the beamformer solution. That is usually done using the **[ft_prepare_leadfield](https://github.com/fieldtrip/fieldtrip/blob/release/ft_prepare_leadfield.m)** function, and the main purpose is to speed up the source reconstruction. The output of **[ft_prepare_leadfield](https://github.com/fieldtrip/fieldtrip/blob/release/ft_prepare_leadfield.m)** can be used as input in **[ft_sourceanalysis](https://github.com/fieldtrip/fieldtrip/blob/release/ft_sourceanalysis.m)**, and you should construct a MATLAB structure that resembles the output of **[ft_prepare_leadfield](https://github.com/fieldtrip/fieldtrip/blob/release/ft_prepare_leadfield.m)**. The best way of getting to know the format of the precomputed leadfields is by computing the leadfield for a spherical model and look at it.
%
% That can be demonstrated using the following code
%
% construct some random EEG electrodes on a sphere
elec     = [];
elec.pnt = randn(128,3);
dum = sqrt(sum(elec.pnt.^2,2));
elec.pnt = 12 * elec.pnt ./ [dum dum dum];  % scale them to a 12 cm sphere
for i=1:128
 elec.label{i} = sprintf('%03d', i);
end

% I suggest that you do not use these electrodes, but instead that
% you get more meaningful electrode locations. See for example http://robertoostenveld.nl/?p=5
%
% create a concentric 3-sphere volume conductor, the radius is the same as for the electrodes
vol   = [];
vol.r = 12 * [0.88 0.92 1.00]; % radii of spheres, the head radius is 12 cm
vol.c = [1 1/80 1];            % conductivity
vol.o = [0 0 0];               % center of sphere

% Using the geometrical description of the sensor locations and the geometrical description of the volume conductor model, you now can compute the leadfield for the locations of interest. By default we scan with the beamformer on a regularly spaced 3D grid of approximately 1 cm resolution. Here is an example for computing the leadfields.
%
% compute the leadfields that can be used for the beamforming
cfg                 = [];
cfg.elec            = elec;
cfg.headmodel       = vol;
cfg.resolution = 2;
cfg.unit       = 'cm';   % same unit as above, i.e. in cm
grid = ft_prepare_leadfield(cfg);

% The idea now is that, using your own forward model, you construct a "grid" structure with the same elements. Note that the channel numbering should correspond with the channels in the data, and that the EEG leadfields should be averaged referenced.
