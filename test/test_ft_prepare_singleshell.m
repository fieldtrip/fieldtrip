function test_ft_prepare_singleshell

% MEM 1gb
% WALLTIME 00:03:29

% TEST test_ft_prepare_singleshell
% TEST ft_headmodel_singleshell ft_prepare_headmodel ft_prepare_singleshell

% function to test ft_headmodel_singleshell. this function is called
% by ft_prepare_headmodel

% the function should work either on an input (segmented) mri, or it
% should have a description of the geometry in the input, or it should
% have a hdmfile (string) that specifies which file to read

success = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the data which is needed

% read in the segmented mri
cd('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
load segmentedmri
mri = segmentedmri; clear segmentedmri;

% specify the file for the headshape
hdmfile  = '/home/common/matlab/fieldtrip/data/Subject01.shape';

% read in the headshape
shape = ft_read_headshape(hdmfile);

%%%%%%%%%%%%%%%%%%%%%
% do the computations

% with an MRI in the input
cfg         = [];
cfg.method  = 'singleshell';
vol1        = ft_prepare_headmodel(cfg, mri);

cfg         = [];
vol1b       = ft_prepare_singleshell(cfg, mri);

% the following needs to be done to be able to make the comparison
vol1 = rmfield(vol1, 'cfg');
vol1b = rmfield(vol1b,'cfg');
vol1bu=ft_convert_units(vol1b,vol1.unit);

% Reason for dashboard failure:
% e.g.  vol1.bnd.pnt(5,:) is different than vol1bu.bnd.pnt(5,:)
% e.g.  vol1.bnd.pnt(8,:) is different than vol1bu.bnd.pnt(8,:)

success     = success && isequal(vol1, vol1bu);
if ~success
  error('ft_prepare_singleshell and ft_prepare_headmodel gave different outputs');
end

% with a filename in the input
cfg         = [];
cfg.method  = 'singleshell';
% cfg.headshape = hdmfile;
% vol2        = ft_prepare_headmodel(cfg);
vol2        = ft_prepare_headmodel(cfg, shape);

cfg         = [];
cfg.headshape = hdmfile;
vol2b       = ft_prepare_singleshell(cfg);

% the following needs to be done to be able to make the comparison
vol2 = rmfield(vol2, 'cfg');
vol2b = rmfield(vol2b,'cfg');

% note from johzum: not sure if vol2 and vol2b should be equal, since the
% old way (vol2b) calls ft_prepare_mesh which creates .tri whereas the new 
% version (vol2) does not have this.  Perhaps only vol2.bnd.pnt should be 
% compared with vol2b.bnd.pnt?
success     = success && isequal(vol2, vol2b);
success     = success && isequal(vol2.bnd.pnt, vol2b.bnd.pnt);
% Note from johzum: the above line also fails. The values are not near each
% other at all, and it's not a units or machine precision issue.
if ~success
  error('ft_prepare_singleshell and ft_prepare_headmodel gave different outputs');
end

% with a point cloud in the input
cfg         = [];
cfg.method  = 'singleshell';
cfg.geom    = shape;
vol3        = ft_prepare_headmodel(cfg);

cfg         = [];
cfg.headshape = shape;
vol3b       = ft_prepare_singleshell(cfg);

% the following needs to be done to be able to make the comparison
vol3 = rmfield(vol3, 'cfg');
vol3b = rmfield(vol3b,'cfg');

success     = success && isequal(vol3, vol3b);
if ~success
  error('ft_prepare_singleshell and ft_prepare_headmodel gave different outputs');
end
