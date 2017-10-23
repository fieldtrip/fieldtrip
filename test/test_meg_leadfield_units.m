function test_meg_leadfield_units

% MEM 2000mb
% WALLTIME 00:10:00

% TEST ft_convert_units ft_datatype_sens ft_convert_vol_sens ft_compute_leadfield current_dipole

%% do a forward computation for a single vector magnetometer
clear all

grad = [];
grad.coilpos = [0 0 12; 0 0 12; 0 0 12];
grad.coilori = eye(3);
grad.label = {'Hx'; 'Hy'; 'Hz'};
grad.unit = 'cm';

% the implementation for the current dipole in an homogenous infinite medium
% returns correct values for inputs that are in SI units
vol0 = [];
vol0.type = 'infinite_currentdipole';
vol0.unit = 'm';

vol1 = [];
vol1.type = 'singlesphere';
vol1.unit = 'cm';
vol1.r = 10;
vol1.o = [0 0 0];

grad = ft_datatype_sens(grad);
grad = ft_convert_units(grad, 'm');
vol0 = ft_convert_units(vol0, 'm');
vol1 = ft_convert_units(vol1, 'm');

dip = [0 0 0.08]; % in meter

lf0 = ft_compute_leadfield(dip, grad, vol0);
lf1 = ft_compute_leadfield(dip, grad, vol1);

n0 = norm(lf0);
n1 = norm(lf1);
% it should be exactly 1/3
assert(n1/n0 < 0.35);
assert(n1/n0 > 0.31);

%% do a forward computation for a more realistic sensor layout
clear all

[pnt, tri] = icosahedron162;
sel = find(pnt(:,3)>0);

grad = [];
grad.pnt = pnt(sel,:) .* 12; % distributed on a 12 cm sphere
grad.ori = pnt(sel,:);
grad.tra = eye(length(sel));
for i=1:length(sel)
  grad.label{i} = sprintf('magnetometer%d', i);
end
grad.unit = 'cm';

vol0 = [];
vol0.type = 'infinite_currentdipole';
vol0.unit = 'cm';

vol1 = [];
vol1.type = 'singlesphere';
vol1.unit = 'cm';
vol1.r = 10;
vol1.o = [0 0 0];

grad = ft_datatype_sens(grad);
grad = ft_convert_units(grad, 'm');
vol0 = ft_convert_units(vol0, 'm');
vol1 = ft_convert_units(vol1, 'm');

dip = [0 0 0.08]; % in meter

% this is to make a selection of the MEG channels
[vol0, grad] = ft_prepare_vol_sens(vol0, grad);
[vol1, grad] = ft_prepare_vol_sens(vol1, grad);

lf0 = ft_compute_leadfield(dip, grad, vol0);
lf1 = ft_compute_leadfield(dip, grad, vol1);

% it should be exactly one
n0 = norm(lf0);
n1 = norm(lf1);
assert(n1/n0 < 1.01);
assert(n1/n0 > 0.99);

figure
ft_plot_dipole(dip, [1 0 0], 'unit', 'm');
ft_plot_vol(vol1);
ft_plot_sens(grad, 'coildiameter', 0.01); % 10 mm
ft_plot_topo3d(grad.chanpos, lf1(:,1));
alpha 0.5

figure
subplot(2,2,1);
plot([lf0(:,1) lf1(:,1)]);
subplot(2,2,2);
plot([lf0(:,2) lf1(:,2)]);
subplot(2,2,3);
plot([lf0(:,3) lf1(:,3)]);


%% do a forward computation for a CTF151 sensor layout
clear all

grad = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds'), 'senstype', 'meg');

vol0 = [];
vol0.type = 'infinite_currentdipole';
vol0.unit = 'cm';

vol1 = [];
vol1.type = 'singlesphere';
vol1.unit = 'cm';
vol1.r = 10;
vol1.o = [0 0 4];

grad = ft_datatype_sens(grad);
grad = ft_convert_units(grad, 'm');
vol0 = ft_convert_units(vol0, 'm');
vol1 = ft_convert_units(vol1, 'm');

dip = [0 0 0.12]; % in meter

% this is to make a selection of the MEG channels
[vol0, grad] = ft_prepare_vol_sens(vol0, grad, 'channel', ft_channelselection('MEG', grad.label));
[vol1, grad] = ft_prepare_vol_sens(vol1, grad, 'channel', ft_channelselection('MEG', grad.label));

lf0 = ft_compute_leadfield(dip, grad, vol0);
lf1 = ft_compute_leadfield(dip, grad, vol1);

n0 = norm(lf0);
n1 = norm(lf1);
% it should be close to one
assert(n1/n0 < 1.1);
assert(n1/n0 > 0.9);

figure
ft_plot_dipole(dip, [1 0 0], 'unit', 'm');
ft_plot_vol(vol1);
ft_plot_sens(grad, 'coildiameter', 0.01); % 10 mm
ft_plot_topo3d(grad.chanpos, lf1(:,1));
alpha 0.5

figure
subplot(2,2,1);
plot([lf0(:,1) lf1(:,1)]);
subplot(2,2,2);
plot([lf0(:,2) lf1(:,2)]);
subplot(2,2,3);
plot([lf0(:,3) lf1(:,3)]);

%% do a forward computation for a CTF151 sensor layout
clear all

grad = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds'), 'senstype', 'meg');

[pnt, tri] = icosahedron162;

mesh = [];
mesh.pnt = pnt * 10;                % set of points on a 10 cm sphere
mesh.pnt(:,3) = mesh.pnt(:,3) + 4;  % shift upward by 4 cm
mesh.unit = 'cm';

vol0 = [];
vol0.type = 'infinite_currentdipole';
vol0.unit = 'cm';

cfg = [];
cfg.method = 'singlesphere';
vol1 = ft_prepare_headmodel(cfg, mesh);

cfg = [];
cfg.method = 'localspheres';
cfg.grad = grad;
vol2 = ft_prepare_headmodel(cfg, mesh);

cfg = [];
cfg.method = 'singleshell';
vol3 = ft_prepare_headmodel(cfg, mesh);

grad = ft_datatype_sens(grad);
grad = ft_convert_units(grad, 'm');
vol0 = ft_convert_units(vol0, 'm');
vol1 = ft_convert_units(vol1, 'm');
vol2 = ft_convert_units(vol2, 'm');
vol3 = ft_convert_units(vol3, 'm');

dip = [0 0 0.12]; % in meter

% this is to make a selection of the MEG channels
[vol0, grad] = ft_prepare_vol_sens(vol0, grad, 'channel', ft_channelselection('MEG', grad.label));
[vol1, grad] = ft_prepare_vol_sens(vol1, grad, 'channel', ft_channelselection('MEG', grad.label));
[vol2, grad] = ft_prepare_vol_sens(vol2, grad, 'channel', ft_channelselection('MEG', grad.label));
[vol3, grad] = ft_prepare_vol_sens(vol3, grad, 'channel', ft_channelselection('MEG', grad.label));

lf0 = ft_compute_leadfield(dip, grad, vol0);
lf1 = ft_compute_leadfield(dip, grad, vol1);
lf2 = ft_compute_leadfield(dip, grad, vol2);
lf3 = ft_compute_leadfield(dip, grad, vol3);

n0 = norm(lf0);
n1 = norm(lf1);
n2 = norm(lf2);
n3 = norm(lf3);
% these should be close to one
assert(abs(n1/n0-1)<0.05);
assert(abs(n2/n0-1)<0.05);
assert(abs(n3/n0-1)<0.05);

figure
ft_plot_dipole(dip, [1 0 0], 'unit', 'm');
ft_plot_vol(vol3);
ft_plot_sens(grad, 'coilsize', 0.01); % 10 mm
ft_plot_topo3d(grad.chanpos, lf3(:,1));
alpha 0.5


