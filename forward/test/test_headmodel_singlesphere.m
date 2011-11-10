function test_headmodel_singlesphere

% TEST test_headmodel_singlesphere
% TEST ft_prepare_vol_sens ft_compute_leadfield ft_prepare_headmodel ft_headmodel_singlesphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first part serves to test basic functionality i.e. is a sphere produced

% function to test ft_headmodel_singlesphere. this function is called by
% ft_prepare_headmodel

% the function should work either on an input (segmented) mri, or it should
% have a description of the geometry in the input, or it should have a
% hdmfile (string) that specifies which file to read

% read in the segmented mri
cd('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer');
load segmentedmri
mri = segmentedmri; clear segmentedmri;

hdmfile = '/home/common/matlab/fieldtrip/data/Subject01.shape';
shape   = ft_read_headshape(hdmfile);

cfg         = [];
cfg.method  = 'singlesphere';
vol1        = ft_prepare_headmodel(cfg, mri);

cfg         = [];
cfg.method  = 'singlesphere';
cfg.hdmfile = hdmfile;
vol2        = ft_prepare_headmodel(cfg);

cfg         = [];
cfg.method  = 'singlesphere';
cfg.geom    = shape;
vol3        = ft_prepare_headmodel(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the second part generates leadfields

addpath('/home/common/matlab/fieldtrip_private'); % in order to find icosahedron

% generate some random points and scale to a unit sphere
pnt = randn(100,3);
for i=1:length(pnt)
  pnt(i,:) = pnt(i,:) / norm(pnt(i,:));
end

geom = [];
geom.pnt = pnt;
vol = ft_headmodel_singlesphere(geom, 'conductivity', [0.42]);

% generate a unit sphere
[pnt, tri] = icosahedron162;

% create a set of electrodes
sel = find(pnt(:,3)>0);
sens.elecpos = pnt(sel,:) * 100;
sens.chanpos = pnt(sel,:) * 100;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% conpute an example leadfield
lf = ft_compute_leadfield([0 0 50], sens, vol);

% figure;
% subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf(:,1))
% subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf(:,2))
% subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf(:,3))

