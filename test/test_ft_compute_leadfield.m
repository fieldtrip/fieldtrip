function test_ft_compute_leadfield

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_compute_leadfield eeg_halfspace_medium_leadfield eeg_halfspace_monopole eeg_leadfieldb eeg_slab_monopole inf_medium_leadfield leadfield_fns leadfield_simbio magnetic_dipole meg_forward meg_leadfield1 ft_prepare_vol_sens
% DATA private

% this function should test all methods

% FIXME eeg

% meg singlesphere
vol = ft_read_headmodel(dccnpath('/project/3031000.02/test/latest/vol/Subject01vol_singlesphere.mat'));
load(dccnpath('/project/3031000.02/test/latest/raw/meg/preproc_ctf151.mat'));
sens = data.grad;
[vol, sens] = ft_prepare_vol_sens(vol, sens);
lf = ft_compute_leadfield([4 3 8], sens, vol);

% meg singleshell
vol = ft_read_headmodel(dccnpath('/project/3031000.02/test/latest/vol/Subject01vol_singleshell.mat'));
load(dccnpath('/project/3031000.02/test/latest/raw/meg/preproc_ctf151.mat'));
sens = data.grad;
[vol, sens] = ft_prepare_vol_sens(vol, sens);
lf = ft_compute_leadfield([4 3 8], sens, vol);

% meg localsphere
vol = ft_read_headmodel(dccnpath('/project/3031000.02/test/latest/vol/Subject01vol_localspheres.mat'));
load(dccnpath('/project/3031000.02/test/latest/raw/meg/preproc_ctf151.mat'));
sens = data.grad;
[vol, sens] = ft_prepare_vol_sens(vol, sens);
lf = ft_compute_leadfield([4 3 8], sens, vol);
