function test_bug1901

% MEM 2500mb
% WALLTIME 00:15:00

% TEST ft_prepare_leadfield ft_prepare_sourcemodel prepare_headmodel ft_convert_units

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% this is some data that should be relatively compatible with the original data
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1901.mat'), 'vol'); 
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1901.mat'), 'grad');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start with the arbitrary units in the vol and grad
cfg = [];
cfg.reducerank = 2;
cfg.feedback = 'off';
cfg.inwardshift = 0;
cfg.grad = grad;
cfg.vol = vol;
cfg.channel = ft_channelselection('MEG', grad.label);
cfg.grid.resolution = 5;
cfg.grid.unit = 'cm';

grid = ft_prepare_leadfield(cfg);

% display the structures
grad;
vol
grid;

% the grid is in cm, which corresponds to the units of the grad, not the vol
% with a 5 cm grid, you can fit 12 sources in the head
assert(sum(grid.inside)==12, 'expected 12 sources inside the volume conductor');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use ft_prepare_sourcemodel instead of ft_prepare_leadfield to speed it up

% repeat it for explicit mm units
vol_mm = ft_convert_units(vol, 'mm');
grad_mm = ft_convert_units(grad, 'mm');

cfg = [];
cfg.reducerank = 2;
cfg.feedback = 'off';
cfg.inwardshift = 0;
cfg.grad = grad_mm;
cfg.vol = vol_mm;
cfg.channel = ft_channelselection('MEG', grad.label);
cfg.grid.resolution = 5; % this is now in mm
grid_mm = ft_prepare_sourcemodel(cfg);
cfg.sourceunits = 'cm'; 
grid_mm2 = ft_prepare_sourcemodel(cfg);

assert(sum(grid_mm.inside)>1e4 && sum(grid_mm.inside)<1e5); % expecting 11822 inside grid points
assert(sum(grid_mm2.inside)==12);                           % expecting 12 inside grid points

% and for explicit cm units
vol_cm = ft_convert_units(vol, 'cm');
grad_cm = ft_convert_units(grad, 'cm');

cfg = [];
cfg.reducerank = 2;
cfg.feedback = 'off';
cfg.inwardshift = 0;
cfg.grad = grad_cm;
cfg.vol = vol_cm;
cfg.channel = ft_channelselection('MEG', grad.label);
cfg.grid.resolution = 5;  % this is now in mm
grid_cm = ft_prepare_sourcemodel(cfg);
cfg.sourceunits = 'mm'; 
grid_cm2 = ft_prepare_sourcemodel(cfg);

assert(sum(grid_cm.inside)==12);                              % expecting 12 inside grid points
assert(sum(grid_cm2.inside)>1e4 && sum(grid_cm2.inside)<1e5); % expecting 11822 inside grid points

% and for explicit m units
vol_m = ft_convert_units(vol, 'm');
grad_m = ft_convert_units(grad, 'm');

cfg = [];
cfg.reducerank = 2;
cfg.feedback = 'off';
cfg.inwardshift = 0;
cfg.grad = grad_m;
cfg.vol = vol_m;
cfg.channel = ft_channelselection('MEG', grad.label);
cfg.grid.resolution = 5; % this is now in m
grid_m = ft_prepare_sourcemodel(cfg);
cfg.sourceunits = 'cm'; 
grid_m2 = ft_prepare_sourcemodel(cfg);

assert(sum(grid_m.inside)==0);   % expecting zero inside grid points
assert(sum(grid_m2.inside)==12); % expecting 12 inside grid points


