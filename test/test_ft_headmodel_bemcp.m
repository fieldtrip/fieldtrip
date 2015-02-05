function test_ft_headmodel_bemcp

% MEM 12gb (or less?)
% WALLTIME 03:00:00 (or less?)

% TEST test_ft_prepare_bemcp
% TEST ft_headmodel_localspheres ft_prepare_localspheres

% to test actual numerical output of bemcp, rather than simply successful
% running and correct inputs (which is what test_ft_prepare_headmodel tests).

% read in the mri
load standard_mri
% this is already (non?)linearly aligned with MNI

% flag for location where running code
atdonders=1;

  %% method 1: segment MRI then compute bnd from it

cfg           = [];
cfg.output    = {'brain','skull','scalp'};
mysegmentedmri  = ft_volumesegment(cfg, mri);

cfg=[];
cfg.tissue={'brain','skull','scalp'};
cfg.numvertices = [1500 1000 500]; % these numbers later match the precomputed bnd
bnd_new=ft_prepare_mesh(cfg,mysegmentedmri);

% I'm choosing 'mm' as units since that is what both the MRI and the
% standard_bem are in already.  Whether that is optimal for 'bemcp' is not
% clear to me.


% Create simple concentricspheres to compare to BEMCP
cfg        = [];
cfg.method ='concentricspheres';
vol_new_cs        = ft_prepare_headmodel(cfg, ft_convert_units(bnd_new,'mm'));
cfg        = [];
cfg.method ='bemcp';
vol_new_bemcp        = ft_prepare_headmodel(cfg, ft_convert_units(bnd_new,'mm'));

%% method 2: compute bnd from segmented MRI

if 0 % I would like to test this, but the MRI is not coregistered and not sure how to coreg it without .anatomy
  
  % read in the already-segmented mri (in case its segmentation is better
  % than above)
  if atdonders
    load('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_eeg/segmentedmri.mat')
  else
    load('/home/zumerj/home/fieldtrip_svn/data/segmentedmri.mat');
  end
  
  % can't do this as no .anatomy
  cfg=[];
  cfg.coordsys='ras';
  cfg.nonlinear='no';
  ft_volumenormalise(cfg,segmentedmri)
  
  
  cfg=[];
  cfg.tissue={'brain','skull','scalp'};
  cfg.numvertices = [3000 2000 1000];
  bnd_seg=ft_prepare_mesh(cfg,segmentedmri);
  
  figure;ft_plot_mesh(bnd_new(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
  hold on;ft_plot_mesh(bnd_seg(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
  
  
  % Create simple concentricspheres to compare to BEMCP
  cfg        = [];
  cfg.method ='concentricspheres';
  vol_seg_cs        = ft_prepare_headmodel(cfg, ft_convert_units(bnd_seg,'mm'));
  cfg        = [];
  cfg.method ='bemcp';
  vol_seg_bemcp        = ft_prepare_headmodel(cfg, ft_convert_units(bnd_seg,'mm'));
end


%% method 3: use existing bnd in standard_bem

% presumably the mesh that was used to create this 'standard' output from
% dipoli should work (irrespective of segmentation above).
load standard_bem
bnd_exist=vol.bnd;
clear vol


figure;ft_plot_mesh(bnd_new(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;ft_plot_mesh(bnd_exist(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
% aside from different order of bnd, they align spatially

% Create simple concentricspheres to compare to BEMCP
cfg        = [];
cfg.method ='concentricspheres';
vol_exist_cs        = ft_prepare_headmodel(cfg, ft_convert_units(bnd_exist,'mm'));
cfg        = [];
cfg.method ='bemcp';
vol_exist_bemcp_mm        = ft_prepare_headmodel(cfg, ft_convert_units(bnd_exist,'mm'));

cfg        = [];
cfg.method ='bemcp';
vol_exist_bemcp_cm        = ft_prepare_headmodel(cfg, ft_convert_units(bnd_exist,'cm'));

cfg        = [];
cfg.method ='bemcp';
vol_exist_bemcp_m        = ft_prepare_headmodel(cfg, ft_convert_units(bnd_exist,'m'));

%% Compare dipoli also

load('standard_bem');
vol_exist_dipoli=vol;
clear vol

%% Compute leadfields

load standard_sourcemodel3d8mm;
sourcemodel=ft_convert_units(sourcemodel,'mm');
if atdonders
  elec=ft_read_sens('/home/common/matlab/fieldtrip/template/electrode/standard_1005.elc');
else
  elec=ft_read_sens('~/home/fieldtrip_svn/template/electrode/standard_1005.elc');
end
elec=ft_convert_units(elec,'mm');

% check for general spatial alignment
figure;ft_plot_mesh(vol_new_bemcp.bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;ft_plot_sens(elec);
hold on;ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));

figure;ft_plot_mesh(vol_exist_bemcp.bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;ft_plot_sens(elec);
hold on;ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));

% make sure vol*.mat of bemcp not NaN
% see bug 1954
if any(isnan(vol_exist_bemcp_mm.mat(:)))
  warning('NaN in vol_exist_bemcp_mm.mat')
end
if any(isnan(vol_new_bemcp.mat(:)))
  warning('NaN in vol_new_bemcp.mat')
end

cfg=[];
cfg.grid=sourcemodel;
cfg.elec=elec;

cfg.vol=vol_new_cs;
lf_new_cs=ft_prepare_leadfield(cfg);

cfg.vol=vol_new_bemcp;
lf_new_bemcp=ft_prepare_leadfield(cfg);

if 0
  cfg.vol=vol_seg_cs;
  lf_seg_cs=ft_prepare_leadfield(cfg);
  
  cfg.vol=vol_seg_cp;
  lf_seg_bemcp=ft_prepare_leadfield(cfg);
end

cfg.vol=vol_exist_cs;
lf_exist_cs=ft_prepare_leadfield(cfg);

cfg.vol=vol_exist_bemcp_mm;
lf_exist_bemcp_mm=ft_prepare_leadfield(cfg);

cfg.vol=vol_exist_dipoli;
lf_exist_dipoli=ft_prepare_leadfield(cfg);

cfg=[];
cfg.vol=vol_exist_bemcp_cm;
cfg.grid=ft_convert_units(sourcemodel,'cm');
cfg.elec=ft_convert_units(elec,'cm');
lf_exist_bemcp_cm=ft_prepare_leadfield(cfg);

cfg=[];
cfg.vol=vol_exist_bemcp_m;
cfg.grid=ft_convert_units(sourcemodel,'m');
cfg.elec=ft_convert_units(elec,'m');
lf_exist_bemcp_m=ft_prepare_leadfield(cfg);


% view results
% Hack: to make LF appear as if it were an ERP

if 0 % only run this manually, not as part of test-battery
  clear lfnames lf
  lfnames=whos('lf*');
  for ll=1:length(lfnames)
    lf=eval(lfnames(ll).name);
    tlock=[];
    tlock.time=[1 2 3];
    tlock.avg=lf.leadfield{dsearchn(sourcemodel.pos,[-20 0 50])};
    tlock.label=lf.cfg.channel;
    tlock.dimord='chan_time';
    
    cfg=[];
    cfg.layout='elec1010.lay';
    cfg.xlim=[0.9 1.1];
    figure;
    ft_topoplotER(cfg,tlock);
    cfg.xlim=[1.9 2.1];
    figure;
    ft_topoplotER(cfg,tlock);
    cfg.xlim=[2.9 3.1];
    figure;
    ft_topoplotER(cfg,tlock);
    
    disp(lfnames(ll).name)
    
    keyboard; % pause to take note of results
    close all;
  end
  % For me, sensible patterns come from the lf*cs and lf*dipoli, but not
  % lf*bemcp*
end
