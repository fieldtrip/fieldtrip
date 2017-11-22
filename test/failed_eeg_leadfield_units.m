function failed_eeg_leadfield_units

% MEM 2000mb
% WALLTIME 00:10:00

% TEST test_meg_leadfield_units
% TEST ft_convert_units ft_datatype_sens ft_convert_vol_sens ft_compute_leadfield current_dipole
[pnt, tri] = icosahedron162;
sel = find(pnt(:,3)>0);

elec = [];
elec.pnt = pnt(sel,:) .* 10; % distributed on a 10 cm sphere
% elec.tra = eye(length(sel));
for i=1:length(sel)
  elec.label{i} = sprintf('electrode%d', i);
end
elec.unit = 'cm';

mesh = [];
mesh(1).pnt = pnt * 10.0; % in cm
mesh(1).tri = tri;
mesh(2).pnt = pnt *  9.5;
mesh(2).tri = tri;
mesh(3).pnt = pnt *  9.0;
mesh(3).tri = tri;

elec = ft_datatype_sens(elec);
elec = ft_convert_units(elec, 'm');
mesh = ft_convert_units(mesh, 'm');

vol0 = [];
vol0.type = 'infinite_currentdipole';
vol0.unit = 'm';

cfg = [];
cfg.method = 'singlesphere';
vol1 = ft_prepare_headmodel(cfg, mesh(1)); % only pass the outermost surface from the mesh

cfg = [];
cfg.method = 'concentricspheres';
cfg.conductivity = [1 1 1]; % functionally identical to a single sphere
vol2 = ft_prepare_headmodel(cfg, mesh);

cfg = [];
cfg.method = 'bemcp';
cfg.conductivity = [1 1 1]; % functionally identical to a single sphere
vol3 = ft_prepare_headmodel(cfg, mesh);

cfg = [];
cfg.method = 'openmeeg';
cfg.conductivity = [1 1 1]; % functionally identical to a single sphere
vol4 = ft_prepare_headmodel(cfg, mesh);

dip = [0 0 0.08]; % in meter

% this is to make a selection of the MEG channels
[vol0, elec] = ft_prepare_vol_sens(vol0, elec);
[vol1, elec] = ft_prepare_vol_sens(vol1, elec);
[vol2, elec] = ft_prepare_vol_sens(vol2, elec);
[vol3, elec] = ft_prepare_vol_sens(vol3, elec);
[vol4, elec] = ft_prepare_vol_sens(vol4, elec);

lf0 = ft_compute_leadfield(dip, elec, vol0);
lf1 = ft_compute_leadfield(dip, elec, vol1);
lf2 = ft_compute_leadfield(dip, elec, vol2);
lf3 = ft_compute_leadfield(dip, elec, vol3);
lf4 = ft_compute_leadfield(dip, elec, vol4);

n0 = norm(lf0);
n1 = norm(lf1);
n2 = norm(lf2);
n3 = norm(lf3);
n4 = norm(lf4);
% these should be relatively close to two, due to the mirror sources for the boundary condition
assert(abs(n1/n0-2)<0.3);
assert(abs(n2/n0-2)<0.3);
assert(abs(n3/n0-2)<0.3);
assert(abs(n3/n0-2)<0.3);
assert(abs(n4/n0-2)<0.3);

figure
ft_plot_dipole(dip, [1 0 0], 'unit', 'm');
ft_plot_vol(vol3, 'edgecolor', 'none', 'facealpha', 0.2);
ft_plot_sens(elec);
ft_plot_topo3d(elec.chanpos, lf3(:,1), 'facealpha', 0.6);
