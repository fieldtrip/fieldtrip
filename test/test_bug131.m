function test_bug131

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_prepare_leadfield

% test the issue related to the scaling of the leadfields in the different implementations

[pnt, tri] = mesh_sphere(162);

% create volume conductor models
vol = [];
vol.o = [0 0 0];
vol.r = 8;
vol.unit = 'm';
vol.type = 'singlesphere';

vol2.bnd.pnt = pnt.*8;
vol2.bnd.tri = tri;
vol2.unit = 'm';
vol2.type = 'singleshell';

% create sensor array
nrm = normals(pnt,tri,'vertex');
grad.pnt = pnt.*10;
grad.pnt(pnt(:,3)<0,:) = [];
grad.ori = nrm(pnt(:,3)>=0,:);
grad.tra = eye(size(grad.pnt,1));
for k = 1:size(grad.pnt,1)
  grad.label{k,1} = ['chan',num2str(k,'%03d')];
end
grad.unit = 'm';

% create dipole grid
grid = [];
sourcemodel.pos = [0 0 4];
sourcemodel.inside = 1;
sourcemodel.outside = [];

% create leadfield with single sphere
cfg = [];
cfg.headmodel = vol;
cfg.sourcemodel = grid;
cfg.grad = grad;
grid1 = ft_prepare_leadfield(cfg);

% create leadfield with singleshell
cfg = [];
cfg.headmodel = vol2;
cfg.sourcemodel = grid;
cfg.grad = grad;
grid2 = ft_prepare_leadfield(cfg);

lf1 = grid1.leadfield{1};
lf2 = grid2.leadfield{1};

% observation: without Gareth's change, the leadfields differ in magnitude
% on the order of 1e10, so the change actually equalizes the magnitudes.
% yet, this in my understanding then only holds for geometrical objects
% defined in SI-units, i.e. in meters. This should then be enforced by the
% higher level function to be able to interpret the units correctly
