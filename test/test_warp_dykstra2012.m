function test_warp_dykstra2012

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY

% function to test the back projection of electrodes using warp_dykstra2012
% Arjen Stolk, Dec 2018

% 'headshape' (random mesh)
x = randn(10,1);
y = randn(10,1);
z = randn(10,1);
shape.tri = delaunay(x, y, z);
shape.pos = [x y z];

%% 2 x 6 grid
elec.label = {'ROFC1';'ROFC2';'ROFC3';'ROFC4';'ROFC5';'ROFC6';'ROFC7';'ROFC8';'ROFC9';'ROFC10';'ROFC11';'ROFC12'};
elec.elecpos = [1 1 1; 2 1 1; 3 1 1; 4 1 1; 5 1 1; 6 1 1; 1 1 2; 2 1 2; 3 1 2; 4 1 2; 5 1 2; 6 1 2];
elec.chanpos = elec.elecpos;
elec.unit = 'mm';

% snap using default method (pos)
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'pos';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end

% snap using label method
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'label';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end

%% 2 x 6 grid with electrode 5 cut out
elec_cutout = elec;
elec_cutout.label(5) = [];
elec_cutout.elecpos(5,:) = [];
elec_cutout.chanpos = elec_cutout.elecpos;

% snap using default method (pos)
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec_cutout;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'pos';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end

% snap using label method
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec_cutout;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'label';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end

%% 2 x 6 grid with electrode 1 cut out
elec_cutout = elec;
elec_cutout.label(1) = [];
elec_cutout.elecpos(1,:) = [];
elec_cutout.chanpos = elec_cutout.elecpos;

% snap using default method (pos)
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec_cutout;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'pos';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end

% snap using label method
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec_cutout;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'label';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end

%% 2 x 6 grid with electrode 12 cut out
elec_cutout = elec;
elec_cutout.label(12) = [];
elec_cutout.elecpos(12,:) = [];
elec_cutout.chanpos = elec_cutout.elecpos;

% snap using default method (pos)
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec_cutout;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'pos';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end

% snap using label method
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec_cutout;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'label';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end

%% 2 x 6 grid with electrodes labeled chan009 and up, misordered (chan011, chan010), and 2 missing (chan014, chan017)
elec_odd.label = {'chan009';'chan011';'chan010';'chan012';'chan013';'chan015';'chan016';'chan018';'chan019';'chan020'};
elec_odd.elecpos = [1 1 1; 2 1 1; 3 1 1; 4 1 1; 5 1 1; 1 1 2; 2 1 2; 4 1 2; 5 1 2; 6 1 2];
elec_odd.chanpos = elec_odd.elecpos;
elec_odd.unit = 'mm';

% snap using default method (pos)
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec_odd;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'pos';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end

% snap using label method
cfg             = [];
cfg.channel     = 'all';
cfg.headshape   = shape;
cfg.elec        = elec_odd;
cfg.method      = 'headshape';
cfg.warp        = 'dykstra2012';
cfg.feedback    = 'yes';
cfg.pairmethod  = 'label';
elec_out = ft_electroderealign(cfg);
if any(any(isnan(elec_out.elecpos)))
  error('the elecpos field should contain non-NaN numbers')
end
