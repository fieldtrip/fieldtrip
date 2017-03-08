function test_headmodel_concentricspheres

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_headmodel_concentricspheres ft_prepare_vol_sens ft_compute_leadfield

% generate some random points
pnt = randn(100,3);
for i=1:length(pnt)
  pnt(i,:) = 100*pnt(i,:) / norm(pnt(i,:));
end
% add some noise
pnt = pnt .* (1 + randn(size(pnt)).*0.1);
geom(1).pnt = pnt;

% generate some random points
pnt = randn(100,3);
for i=1:length(pnt)
  pnt(i,:) = 90*pnt(i,:) / norm(pnt(i,:));
end
% add some noise
pnt = pnt .* (1 + randn(size(pnt)).*0.1);
geom(2).pnt = pnt;


% generate some random points
pnt = randn(100,3);
for i=1:length(pnt)
  pnt(i,:) = 80*pnt(i,:) / norm(pnt(i,:));
end
% add some noise
pnt = pnt .* (1 + randn(size(pnt)).*0.1);
geom(3).pnt = pnt;


% generate some random points
pnt = randn(100,3);
for i=1:length(pnt)
  pnt(i,:) = 75*pnt(i,:) / norm(pnt(i,:));
end
% add some noise
pnt = pnt .* (1 + randn(size(pnt)).*0.1);
geom(4).pnt = pnt;

% the geometry is sorted ourtward in
vol = ft_headmodel_concentricspheres(geom, 'conductivity', [0.42 0.0033 1 0.42]);

% generate a unit sphere
[pnt, tri] = icosahedron162;

% create a set of electrodes
sel = find(pnt(:,3)>0);
sens.elecpos = pnt(sel,:) * 100;
sens.chanpos = pnt(sel,:) * 100;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% conpute an example leadfield
lf = ft_compute_leadfield([0 0 50], sens, vol);

figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf(:,3))
