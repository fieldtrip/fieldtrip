function test_ft_prepare_singleshell

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_headmodel_singleshell ft_prepare_headmodel ft_prepare_singleshell

% function to test ft_headmodel_singleshell. this function is called
% by ft_prepare_headmodel

% the function should work either on an input (segmented) mri, or it
% should have a description of the geometry in the input, or it should
% have a hdmfile (string) that specifies which file to read

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the data which is needed

% read in the segmented mri
datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
load(fullfile(datadir, 'segmentedmri.mat'));
mri = segmentedmri; clear segmentedmri;

% specify the file for the headshape
hdmfile  = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.shape');

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
vol1.bnd    = rmfield(vol1.bnd, 'cfg');
vol1bu.bnd  = rmfield(vol1bu.bnd, 'cfg');

if ~isequal(vol1, vol1bu)
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
vol2  = rmfield(vol2, 'cfg');
vol2b = rmfield(vol2b,'cfg');

vol2.bnd  = rmfield(vol2.bnd, 'cfg');
vol2b.bnd = rmfield(vol2b.bnd,'cfg');

[ok, msg] = isalmostequal(vol2, vol2b, 'abstol', 1e-5, 'diffabs', 1);
  
if ~ok
  disp(msg);
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
vol3  = rmfield(vol3, 'cfg');
vol3b = rmfield(vol3b,'cfg');

vol3.bnd  = rmfield(vol3.bnd, 'cfg');
vol3b.bnd = rmfield(vol3b.bnd,'cfg');

[ok, msg] = isalmostequal(vol3, vol3b, 'abstol', 1e-5, 'diffabs', 1);
if ~ok
  disp(msg);
  error('ft_prepare_singleshell and ft_prepare_headmodel gave different outputs');
end
