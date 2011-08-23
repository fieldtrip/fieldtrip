% script to test ft_headmodel_singlesphere. this function is called by
% ft_prepare_headmodel

% the function should work either on an input (segmented) mri, or it should
% have a description of the geometry in the input, or it should have a
% hdmfile (string) that specifies which file to read

curdir = pwd;

% read in the segmented mri
cd('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
%cd('/Volumes/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
load segmentedmri
mri = segmentedmri; clear segmentedmri;

cfg         = [];
cfg.method  = 'singlesphere';
vol1        = ft_prepare_headmodel(cfg, mri);

cd('/home/common/matlab/fieldtrip/data/');
%cd('/Volumes/home/common/matlab/fieldtrip/data/');
hdmfile     = 'Subject01.shape';

cfg         = [];
cfg.method  = 'singlesphere';
cfg.hdmfile = hdmfile;
vol2        = ft_prepare_headmodel(cfg);

% read in the headshape
shape = ft_read_headshape(cfg.hdmfile);

cfg         = [];
cfg.method  = 'singlesphere';
cfg.geom    = shape;
vol3        = ft_prepare_headmodel(cfg);

cd(curdir);


% % generate some random points and scale to a unit sphere
% pnt = randn(100,3);
% for i=1:length(pnt)
%   pnt(i,:) = pnt(i,:) / norm(pnt(i,:));
% end
% 
% vol = ft_headmodel_singlesphere(pnt, 'conductivity', [0.42]);
% 
% % generate a unit sphere
% [pnt, tri] = icosahedron162;
% 
% % create a set of electrodes
% sel = find(pnt(:,3)>0);
% sens.pnt = pnt(sel,:) * 100;
% for i=1:length(sel)
%   sens.label{i} = sprintf('chan%03d', i);
% end
% 
% % project the electrodes on the volume conduction model
% [vol, sens] = ft_prepare_vol_sens(vol, sens);
% 
% % conpute an example leadfield
% lf = ft_compute_leadfield([0 0 50], sens, vol);
% 
% figure;
% subplot(2,2,1); ft_plot_topo3d(sens.pnt, lf(:,1))
% subplot(2,2,2); ft_plot_topo3d(sens.pnt, lf(:,2))
% subplot(2,2,3); ft_plot_topo3d(sens.pnt, lf(:,3))
% 
