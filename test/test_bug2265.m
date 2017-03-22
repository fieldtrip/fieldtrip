function test_bug2265

% MEM 2000mb
% WALLTIME 00:10:00

% TEST ft_convert_units ft_prepare_sourcemodel

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

clear all
close all

load(dccnpath('/home/common/matlab/fieldtrip/template/headmodel/standard_bem.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/template/headmodel/standard_mri.mat'));
elecs = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/template/electrode/standard_1020.elc'));

%% 
CMMM = 'cm'; scale=1;
CMMM = 'mm'; scale=10; % prepare_leadfield, or another method does not work with mm (a bugg)

vol = ft_convert_units(vol, CMMM); % Convert the vol to cm, or mm

% construct the dipole grid in the template brain coordinates
% the source units are in cm
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg = [];
cfg.grid.xgrid = scale*[-20:1:20];
cfg.grid.ygrid = scale*[-20:1:20];
cfg.grid.zgrid = scale*[-20:1:20];
cfg.grid.unit = CMMM;
cfg.grid.tight = 'yes';
cfg.inwardshift = -scale*0;
cfg.vol = vol;
template_grid = ft_prepare_sourcemodel(cfg);

assert(strcmp(vol.unit, CMMM));
assert(strcmp(template_grid.unit, CMMM));


%% make a figure with the template head model and dipole grid
figure
hold on
ft_plot_vol(vol, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5;
camlight;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

elecs = ft_convert_units(elecs, CMMM);

[vol, elecs] = ft_prepare_vol_sens(vol, elecs);

ft_plot_mesh(elecs.elecpos,'vertexcolor',[1 .3 .3]);


data = [];
cfg.grid = template_grid;
cfg.vol = vol;
cfg.elec = elecs;
cfg.unit = CMMM; % this is confusing in case it is different from cfg.grid.unit
[grid] = ft_prepare_leadfield(cfg, data);

grid = ft_convert_units(grid, CMMM); %% does not really functional!!


%% Dipole simulation
% note that beamformer scanning will be done with a 1cm grid, so you should
% not put the dipole on a position that will not be covered by a grid
% location later
cfg = [];
cfg.vol = vol;
cfg.elec = elecs;
cfg.dip.pos = scale*[
  3 -3 6 % dipole 1
  ];
cfg.dip.mom = [ 1 0 0 ]';

cfg.dip.frequency = [10]*1;
cfg.dip.phase=pi/4*[1];
cfg.dip.amplitude=[10];
cfg.relnoise = .1;
cfg.ntrials = 20;
data = ft_dipolesimulation(cfg);


%% compute the data covariance matrix, which will capture the activity of
% the simulated dipole
cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [data.time{1}(1) data.time{1}(end)];
timelock = ft_timelockanalysis(cfg, data);


%% source analysis
cfg = [];
METHOD = 'lcmv'; % 'lcmv' 'rv' 'mne'
cfg.method = METHOD;
cfg.grid = grid;


cfg.vol = vol;
cfg.elec = elecs;
cfg.snr = 10;
source = ft_sourceanalysis(cfg, timelock);

%% interpolate onto the anatomical MRI

cfg = [];
cfg.downsample = 2;
if strcmp(METHOD,'lcmv')
  cfg.parameter = 'avg.pow';
elseif strcmp(METHOD,'rv')
  cfg.parameter = 'avg.rv';
else
  cfg.parameter = 'avg.pow';
end
source = ft_sourceinterpolate(cfg, source , mri); % Note: this function crashes if mne method is used

%% plot the source reconstruction

cfg = [];
if strcmp(METHOD,'lcmv')
  try, source=rmfield(source,'time'); end
  cfg.funparameter = 'avg.pow';
elseif strcmp(METHOD,'rv')
  try, source=rmfield(source,'time'); end
  cfg.funparameter = 'avg.rv';
else
  cfg.funparameter = 'avg.pow';
end

cfg.method = 'ortho';
ft_sourceplot(cfg, source);
