function test_ft_electroderealign

% MEM 2gb
% WALLTIME 00:03:27

% TEST ft_electroderealign
% TEST ft_read_mri ft_read_sens ft_prepare_mesh warp_apply

% first version of this script was made on 10.03.2013 by Lilla Magyari

interactive = false;  % use this for running it without user interaction 
% interactive = true;  % use this for running it with user interaction 

%% load mri, segmentation and electrode definition
if ispc
    datadir = 'H:';
else
    datadir = '/home';
end
mri=ft_read_mri(strcat(datadir,'/common/matlab/fieldtrip/data/Subject01.mri'));
load(strcat(datadir,'/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_eeg/segmentedmri'));
elec = ft_read_sens(strcat(datadir,'/common/matlab/fieldtrip/template/electrode/standard_1020.elc'));
temp = ft_read_sens(strcat(datadir,'/common/matlab/fieldtrip/template/electrode/standard_1005.elc'));

% create a bem and a fem mesh

cfg=[];
cfg.tissue={'brain','skull','scalp'};
cfg.numvertices=[300 200 100];
bem=ft_prepare_mesh(cfg,segmentedmri);

cfg=[];
cfg.method = 'hexahedral';
cfg.tissue={'brain','skull','scalp'};
fem=ft_prepare_mesh(cfg,segmentedmri);


%% method: fiducial

 nas=mri.hdr.fiducial.mri.nas;
 lpa=mri.hdr.fiducial.mri.lpa;
 rpa=mri.hdr.fiducial.mri.rpa;
 
 transm=mri.transform;
 
 nas=warp_apply(transm,nas, 'homogenous');
 lpa=warp_apply(transm,lpa, 'homogenous');
 rpa=warp_apply(transm,rpa, 'homogenous');
 
 fiducials.chanpos=[nas; lpa; rpa];
 fiducials.label={'Nz','LPA','RPA'};
 fiducials.unit='mm';
 
 
 cfg=[];
 cfg.method='fiducial';
 cfg.template=fiducials;
 cfg.elec = elec;
 cfg.channel = {'Nz', 'LPA', 'RPA','M1','M2'};
 cfg.fiducial={'Nz', 'LPA', 'RPA'};
 elec_align=ft_electroderealign(cfg);
 
%% method: interactive 

if interactive
  cfg=[];
  cfg.method='interactive';
  cfg.elec=elec_align;
  cfg.headshape=bem(3);
%  
  elec_align2=ft_electroderealign(cfg);
  
  cfg=[];
  cfg.method='interactive';
  cfg.elec=elec_align;
  cfg.headshape=fem;
%  
  elec_align3=ft_electroderealign(cfg);
 
 elec_align2=rmfield(elec_align2,'cfg');
 elec_align3=rmfield(elec_align3,'cfg');
  
 assert(isequal(elec_align2,elec_align3),'Oops.');
  
%  close all;
end   
%% method: template - electrode set

cfg=[];
cfg.elec = elec_align;
cfg.template = temp;
cfg.method = 'template';
elec_align4 = ft_electroderealign(cfg);

assert(~isequal(elec_align.chanpos,elec_align4.chanpos));

%% method: template - headshape

% doesn't work for me
% bem
%   cfg=[];
%   cfg.method='template';
%   cfg.elec=elec_align;
%   cfg.headshape=bem(3);
%   elec_align5=ft_electroderealign(cfg);
%   
% % fem
% 
%   cfg=[];
%   cfg.method='template';
%   cfg.elec=elec_align;
%   cfg.headshape=fem;
%   elec_align6=ft_electroderealign(cfg);

%% manual
% doesn't work for me

%   cfg=[];
%   cfg.method='manual';
%   cfg.headshape=bem(3); 
%   elec_align7=ft_electroderealign(cfg);
  

% 
%   cfg=[];
%   cfg.method='manual';
%   cfg.headshape=fem; 
%   elec_align8=ft_electroderealign(cfg);

clear all;
close all;
   
