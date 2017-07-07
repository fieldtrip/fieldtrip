function test_fileformat_asa

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_sens ft_read_vol ft_read_headshape ft_read_mri read_asa read_asa_dip	read_asa_mri read_asa_vol read_asa_bnd	read_asa_elc	read_asa_msr

elcfile  = dccnpath('/home/common/matlab/fieldtrip/data/test/original/electrodes/asa/standard_primed.elc');
volfile  = dccnpath('/home/common/matlab/fieldtrip/data/test/original/headmodel/asa/standard.vol');
skinfile = dccnpath('/home/common/matlab/fieldtrip/data/test/original/headshape/asa/standard_skin_14038.vol');
mrifile  = dccnpath('/home/common/matlab/fieldtrip/data/test/original/mri/asa/standard.mri');

elec = ft_read_sens(elcfile);
vol  = ft_read_vol(volfile);
skin = ft_read_headshape(skinfile);
mri  = ft_read_mri(mrifile);
