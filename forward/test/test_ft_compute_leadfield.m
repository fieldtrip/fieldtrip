function test_ft_compute_leadfield

% TEST test_ft_compute_leadfield
% ft_compute_leadfield eeg_halfspace_medium_leadfield eeg_halfspace_monopole eeg_leadfieldb eeg_slab_monopole inf_medium_leadfield leadfield_fns leadfield_simbio magnetic_dipole meg_forward meg_leadfield1 ft_prepare_vol_sens

% this function should test all methods

% eeg


% meg singlesphere
vol = ft_read_vol('/home/common/matlab/fieldtrip/data/test/vol/Subject01vol_singlesphere.mat');
load('/home/common/matlab/fieldtrip/data/test/raw/meg/preproc_ctf151');
sens = data.grad;
[vol, sens] = ft_prepare_vol_sens(vol, sens);
lf = ft_compute_leadfield([4 3 8], sens, vol);

% meg singleshell
vol = ft_read_vol('/home/common/matlab/fieldtrip/data/test/vol/Subject01vol_singleshell.mat');
load('/home/common/matlab/fieldtrip/data/test/raw/meg/preproc_ctf151');
sens = data.grad;
[vol, sens] = ft_prepare_vol_sens(vol, sens);
lf = ft_compute_leadfield([4 3 8], sens, vol);

% meg multisphere
vol = ft_read_vol('/home/common/matlab/fieldtrip/data/test/vol/Subject01vol_localspheres.mat');
load('/home/common/matlab/fieldtrip/data/test/raw/meg/preproc_ctf151');
sens = data.grad;
[vol, sens] = ft_prepare_vol_sens(vol, sens);
lf = ft_compute_leadfield([4 3 8], sens, vol);
