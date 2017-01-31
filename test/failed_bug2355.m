function failed_bug2355

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug2355
% TEST ft_prepare_sourcemodel

global ft_default
ft_default.checkconfig = 'off';
ft_default.trackconfig = 'off';

%% make a spherical model fitted to scalp surface, in cm
vol = [];
vol.r = 12;
vol.o = [0 0 4];
vol.unit = 'cm';

%% make a set of electrodes, in cm
[pnt, tri] = icosahedron162;
pnt = pnt(pnt(:,3)>0,:); % only above z=0
pnt = pnt*12;
pnt(:,3) = pnt(:,3) + 4;

elec = [];
elec.elecpos = pnt;
for i=1:length(pnt)
  elec.label{i} = num2str(i);
end

%% make 3-D grids as source model
% cfg.sourceunits was used prior to 31 October 2013

clear cfg* grid*

cfg1.elec = elec;
cfg1.vol  = vol;
cfg1.grid.resolution = 1;
cfg1.sourceunits = 'cm';
grid1 = ft_prepare_sourcemodel(cfg1);
assert(sum(grid1.inside)==2469);

cfg2.elec = elec;
cfg2.vol  = ft_convert_units(vol, 'mm');
cfg2.grid.resolution = 1;
cfg2.sourceunits = 'cm';
grid2 = ft_prepare_sourcemodel(cfg2);
assert(sum(grid2.inside)==2469);

cfg3.elec = ft_convert_units(elec, 'mm');
cfg3.vol  = vol;
cfg3.grid.resolution = 1;
cfg3.sourceunits = 'cm';
grid3 = ft_prepare_sourcemodel(cfg3);
assert(sum(grid3.inside)==2469);

cfg4.elec = ft_convert_units(elec, 'mm');
cfg4.vol  = ft_convert_units(vol, 'mm');
cfg4.grid.resolution = 1;
cfg4.sourceunits = 'cm';
grid4 = ft_prepare_sourcemodel(cfg4);
assert(sum(grid4.inside)==2469);

cfg5.elec = ft_convert_units(elec, 'dm');
cfg5.vol  = ft_convert_units(vol, 'mm');
cfg5.grid.resolution = 1;
cfg5.sourceunits = 'cm';
grid5 = ft_prepare_sourcemodel(cfg5);
assert(sum(grid5.inside)==2469);

%% repeat with cfg.xgrid etc. instead of cfg.resolution
% cfg.sourceunits was used prior to 31 October 2013, now it is cfg.grid.unit

clear cfg* grid*

cfg1.elec = ft_convert_units(elec, 'dm');
cfg1.vol  = ft_convert_units(vol, 'mm');
cfg1.resolution = 1;
cfg1.sourceunits = 'cm';
cfg.grid.tight = 'no';
grid1 = ft_prepare_sourcemodel(cfg1);
assert(sum(grid1.inside)==2469);

cfg2.elec = ft_convert_units(elec, 'cm'); % the units will be copied from the sens
cfg2.vol  = ft_convert_units(vol, 'mm');
cfg2.xgrid = -20:1:20;
cfg2.ygrid = -20:1:20;
cfg2.zgrid = 7:1:20; % note that the center of the sphere is not at [0 0 0]
grid2 = ft_prepare_sourcemodel(cfg2);
assert(sum(grid2.inside)==2469);

cfg3.elec = ft_convert_units(elec, 'dm');
cfg3.vol  = ft_convert_units(vol, 'mm');
cfg3.xgrid = -20:1:20;
cfg3.ygrid = -20:1:20;
cfg3.zgrid = 7:1:20; % note that the center of the sphere is not at [0 0 0]
cfg3.sourceunits = 'cm';
grid3 = ft_prepare_sourcemodel(cfg3);
assert(sum(grid3.inside)==2469);

%% repeat with new name of cfg.sourceunits option
% cfg.sourceunits was used prior to 31 October 2013, now it is cfg.grid.unit

clear cfg* grid*

cfg1.elec = elec;
cfg1.vol  = vol;
cfg1.grid.resolution = 1;
cfg1.grid.unit = 'cm';
grid1 = ft_prepare_sourcemodel(cfg1);
assert(sum(grid1.inside)==2469);

cfg2.elec = elec;
cfg2.vol  = ft_convert_units(vol, 'mm');
cfg2.grid.resolution = 1;
cfg2.grid.unit = 'cm';
grid2 = ft_prepare_sourcemodel(cfg2);
assert(sum(grid2.inside)==2469);

cfg3.elec = ft_convert_units(elec, 'mm');
cfg3.vol  = vol;
cfg3.grid.resolution = 1;
cfg3.grid.unit = 'cm';
grid3 = ft_prepare_sourcemodel(cfg3);
assert(sum(grid3.inside)==2469);

cfg4.elec = ft_convert_units(elec, 'mm');
cfg4.vol  = ft_convert_units(vol, 'mm');
cfg4.grid.resolution = 1;
cfg4.grid.unit = 'cm';
grid4 = ft_prepare_sourcemodel(cfg4);
assert(sum(grid4.inside)==2469);

