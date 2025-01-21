function test_fileformat_asa

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_sens ft_read_headmodel ft_read_headshape ft_read_mri read_asa read_asa_dip	read_asa_mri read_asa_vol read_asa_bnd read_asa_elc read_asa_msr
% DATA private

elcfile  = dccnpath('/project/3031000.02/test/original/electrodes/asa/standard_primed.elc');
volfile  = dccnpath('/project/3031000.02/test/original/headmodel/asa/standard.vol');
skinfile = dccnpath('/project/3031000.02/test/original/headshape/asa/standard_skin_14038.vol');
mrifile  = dccnpath('/project/3031000.02/test/original/mri/asa/standard.mri');

elec = ft_read_sens(elcfile);
vol  = ft_read_headmodel(volfile);
skin = ft_read_headshape(skinfile);
mri  = ft_read_mri(mrifile);
