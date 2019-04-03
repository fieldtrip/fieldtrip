% this is the radius and conductivity of each of the three compartments
r = [100 92 88];
c = [1 1/80 1];

% this is the initial description of the spherical mesh
[pnt, tri] = mesh_sphere(162);

% create a set of electrodes on the outer surface
sens.pnt = max(r) * pnt;
sens.label = {};
for i=1:size(pnt,1)
  sens.label{i} = sprintf('vertex%03d', i);
end

% this is the position of the probe dipole
pos = [0 0 70];

% create a BEM volume conduction model
vol = [];
vol.bnd(1).pnt = pnt * r(1);
vol.bnd(1).tri = tri;
vol.bnd(2).pnt = pnt * r(2);
vol.bnd(2).tri = tri;
vol.bnd(3).pnt = pnt * r(3);
vol.bnd(3).tri = tri;
vol.cond = c;
% compute the BEM system matrix
vol = dipoli(vol, true);
[vol, sens] = prepare_vol_sens(vol, sens);

% compute the forward solution
lf = compute_leadfield(pos, sens, vol);

figure; triplot(pnt, tri, lf(:,1));
figure; triplot(pnt, tri, lf(:,2));
figure; triplot(pnt, tri, lf(:,3));

% now show the same example using three concentric spheres
vol_sphere.r = r;
vol_sphere.c = c;

% compute the forward solution
lf_sphere = compute_leadfield(pos, sens, vol_sphere);

figure
plot(lf(:,1))
hold on
plot(lf_sphere(:,1), 'r')
