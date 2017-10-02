function test_bug2338

% MEM 2000mb
% WALLTIME 00:20:00

% TEST ft_prepare_bemmodel ft_prepare_headmodel ft_prepare_leadfield ft_compute_leadfield ft_headmodel_openmeeg 

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% 4 Layers
r = [85 88 92 100];
c = [1 1/20 1/80 1];

order = [1 2 3 4];

% Description of the spherical mesh
[pnt, tri] = icosahedron42;

% Create a set of electrodes on the outer surface
elec.elecpos = max(r) * pnt;
elec.label = {};
nelec = size(elec.elecpos,1);
for ii=1:nelec
  elec.label{ii} = sprintf('vertex%03d', ii);
end

% Create one triangulated mesh for each boundary, the first boundary is inside
mesh = [];
for ii=1:length(r)
  mesh.bnd(ii).pnt = pnt * r(ii);
  mesh.bnd(ii).tri = tri;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the BEM model using the old code
cfg = [];
cfg.method = 'openmeeg';

cfg.conductivity = c(order);
mesh1 = mesh;
mesh1.bnd = mesh1.bnd(order);
vol1 = ft_prepare_bemmodel(cfg, mesh1);

% flip the geometry and conductivity around
order = [4 3 2 1];

cfg.conductivity = c(order);
mesh2 = mesh;
mesh2.bnd = mesh2.bnd(order);
vol2 = ft_prepare_bemmodel(cfg, mesh2);

cfg = [];
cfg.grid.pos = [0 0 70];
cfg.elec = elec;
cfg.vol = vol1;
lf1 = ft_prepare_leadfield(cfg);
cfg.vol = vol2;
lf2 = ft_prepare_leadfield(cfg);

assert(isalmostequal(lf1.leadfield{1}, lf2.leadfield{1}, 'reltol', 1e-6));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the BEM model using the new code
cfg = [];
cfg.method = 'openmeeg';

cfg.conductivity = c(order);
mesh1 = mesh;
mesh1.bnd = mesh1.bnd(order);
vol1 = ft_prepare_headmodel(cfg, mesh1);

% flip the geometry and conductivity around
order = [4 3 2 1];

cfg.conductivity = c(order);
mesh2 = mesh;
mesh2.bnd = mesh2.bnd(order);
vol2 = ft_prepare_headmodel(cfg, mesh2);

cfg = [];
cfg.grid.pos = [0 0 70];
cfg.elec = elec;
cfg.vol = vol1;
lf1 = ft_prepare_leadfield(cfg);
cfg.vol = vol2;
lf2 = ft_prepare_leadfield(cfg);

assert(isalmostequal(lf1.leadfield{1}, lf2.leadfield{1}, 'reltol', 1e-6));

