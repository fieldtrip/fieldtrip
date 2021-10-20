function test_ft_dipolesimulation

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_dipolesimulation

nchan = 32;
% construct a set of electrodes randomly distributed over the upper hemisphere
elec.label = cellstr(num2str((1:nchan).'));
elec.unit = 'cm';
elec.tra = eye(nchan);
elec.elecpos = randn(nchan,3);
elec.elecpos(:,3) = abs(elec.elecpos(:,3));
for i=1:nchan
  elec.elecpos(i,:) = 10*elec.elecpos(i,:)/norm(elec.elecpos(i,:));
end
elec.chanpos = elec.elecpos;

geometry = [];
geometry.pos = elec.elecpos;
geometry.unit = elec.unit;
cfg = [];
cfg.conductivity = [0.33 1.00 0.042 0.33];
headmodel = ft_headmodel_concentricspheres(geometry, 'conductivity', cfg.conductivity);

cfg = [];
cfg.dip.pos = randn(2,3);
cfg.dip.mom = randn(3,2);
cfg.dip.frequency = [10 20];
cfg.dip.phase = [0 45];
cfg.dip.amplitude = [3 5];
cfg.headmodel = headmodel;
cfg.elec = elec;
cfgout = ft_dipolesimulation(cfg);