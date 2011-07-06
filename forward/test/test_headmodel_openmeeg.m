% generate a unit sphere
[pnt, tri] = icosahedron162;

% create the BEM geometries
geom = [];
geom.bnd(1).pnt = pnt * 100;
geom.bnd(1).tri = tri;
geom.bnd(2).pnt = pnt * 90;
geom.bnd(2).tri = tri;
geom.bnd(3).pnt = pnt * 80;
geom.bnd(3).tri = tri;

vol = ft_headmodel_openmeeg(geom, 'conductivity', [1 2 3]);

% create a set of electrodes
sel = find(pnt(:,3)>0);
sens.pnt = pnt(sel,:) * 100;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% conpute an example leadfield
lf = ft_compute_leadfield([0 0 50], sens, vol);

figure;
subplot(2,2,1); ft_plot_topo3d(sens.pnt, lf(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.pnt, lf(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.pnt, lf(:,3))

