function failed_tutorial_headmodel_eeg

% MEM 5gb
% WALLTIME 01:30:00

% TEST test_tutorial_headmodel_eeg
% TEST ft_read_mri ft_volumesegment ft_prepare_mesh ft_prepare_headmodel
% TEST ft_read_sens ft_warp_apply ft_electroderealign
% TEST ft_plot_mesh ft_plot_vol ft_plot_sens

clear all;
%% load mri
mri=ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));

%% segmentation

cfg=[];
cfg.output={'brain','skull','scalp'};
segmentedmri = ft_volumesegment(cfg,mri);

%save segmentedmri segmentedmri;

% check if segmentation is the same as the segmentation on the ftp site

segmentedmri2=load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_eeg/segmentedmri'));
segmentedmri2=rmfield(segmentedmri2.segmentedmri,'cfg');
segmentedmri1=rmfield(segmentedmri,'cfg');
assert(isequal(segmentedmri2,segmentedmri1),'Segmentation does not match the segmentation on ftp/tutorial/headmodel_eeg');
clear segmentedmri1 segmentedmri2;

%% triangulation
%load segmentedmri;

cfg=[];
cfg.tissue={'brain','skull','scalp'};
cfg.numvertices=[3000 2000 1000];
%cfg.sourceunits=segmentedmri.unit;
bnd=ft_prepare_mesh(cfg,segmentedmri);
%save bnd bnd;

clear segmentedmri;
%% headmodel
% dipoli
cfg=[];
cfg.method='dipoli';
vol=ft_prepare_headmodel(cfg,bnd);
%save vol vol;

% Openmeeg
system('module load openmeeg'); % Load openmeeg paths
cfg=[];
cfg.method='openmeeg';
vol_openmeeg=ft_prepare_headmodel(cfg,bnd);
% warning- takes ~40 minutes

% check if segmentation is the same as the segmentation on the ftp site
% dipoli
vol2=load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_eeg/vol'));
vol2=rmfield(vol2.vol,'cfg');
vol1=rmfield(vol,'cfg');
assert(isequal(vol2,vol1),'Segmentation of dipoli vol does not match the segmentation on ftp/tutorial/headmodel_eeg');
clear vol1 vol2;

% openmeeg
vol2=load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_eeg/vol_openmeeg'));
vol2=rmfield(vol2.vol,'cfg');
vol1=rmfield(vol_openmeeg,'cfg');
assert(isequal(vol2,vol1),'Segmentation of openmeeg vol does not match the segmentation on ftp/tutorial/headmodel_eeg');
clear vol1 vol2;


%% visualization
figure;
ft_plot_mesh(vol.bnd(1),'facecolor','none');  %skin
figure;
ft_plot_mesh(vol.bnd(2),'facecolor','none');
figure;
ft_plot_mesh(vol.bnd(3),'facecolor','none');
% 
figure;
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(vol.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]);

close all;
%% align electrodes

elec = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/template/electrode/standard_1020.elc'));

figure;
ft_plot_sens(elec,'style','sk');
hold on;
ft_plot_mesh(vol.bnd(1),'facealpha', 0.85, 'edgecolor', 'none', 'facecolor', [0.65 0.65 0.65]); %scalp
close all;

 nas=mri.hdr.fiducial.mri.nas;
 lpa=mri.hdr.fiducial.mri.lpa;
 rpa=mri.hdr.fiducial.mri.rpa;
 
 transm=mri.transform;
 
 nas=ft_warp_apply(transm,nas, 'homogenous');
 lpa=ft_warp_apply(transm,lpa, 'homogenous');
 rpa=ft_warp_apply(transm,rpa, 'homogenous');
 
 fiducials.chanpos=[nas; lpa; rpa];
 fiducials.label={'Nz','LPA','RPA'};
 fiducials.unit='mm';
 
 % ensure the elec to have a coordsys field, in order to avoid
 % the interactive step due to the call to ft_determine_coordsys
 % in ft_electroderealign
 elec.coordsys = 'ctf';
 
 cfg=[];
 cfg.method='fiducial';
 cfg.template=fiducials;
 cfg.elec = elec;
 cfg.fiducial={'Nz', 'LPA', 'RPA'};
 elec_align=ft_electroderealign(cfg);
 
figure;
ft_plot_sens(elec_align,'style','sk');
hold on;
ft_plot_mesh(vol.bnd(1),'facealpha', 0.85, 'edgecolor', 'none', 'facecolor', [0.65 0.65 0.65]); %scalp
close all;

%   cfg=[];
%   cfg.method='interactive';
%   cfg.elec=elec_align;
%   cfg.headshape=vol.bnd(1);
% %  
%   elec_align=ft_electroderealign(cfg);
%  close all;

  
   
  figure;
ft_plot_sens(elec_align,'style','sk');
hold on;
ft_plot_mesh(vol.bnd(1),'facealpha', 0.85, 'edgecolor', 'none', 'facecolor', [0.65 0.65 0.65]); %scalp

  close all;
%% concentric spheres
% load bnd;
cfg=[];
cfg.method='concentricspheres';
cfg.conductivity=[0.33 0.041 0.33];
vol_cs=ft_prepare_headmodel(cfg,bnd);

figure;
ft_plot_vol(vol_cs,'facealpha', 0.3)
hold on;
ft_plot_mesh(bnd(3),'facealpha', 0.3, 'facecolor', 'red', 'edgecolor', 'none');

figure;
ft_plot_vol(vol_cs,'facealpha', 0.3)
hold on;
ft_plot_mesh(bnd(1),'facealpha', 0.3, 'facecolor', 'red', 'edgecolor', 'none');

close all;
clear all;
