function test_headmodel_singlesphere

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_vol_sens ft_compute_leadfield ft_headmodel_singlesphere

% generate some random points and scale to a unit sphere
pnt = randn(100,3);
for i=1:length(pnt)
  pnt(i,:) = pnt(i,:) / norm(pnt(i,:));
end

geom = [];
geom.pnt = pnt * 100;
vol = ft_headmodel_singlesphere(geom, 'conductivity', [0.42]);

% generate a unit sphere
[pnt, tri] = icosahedron162;

% create a set of electrodes
sel = find(pnt(:,3)>0);
sens.elecpos = pnt(sel,:) * 100;
sens.chanpos = pnt(sel,:) * 100;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end

sens = ft_convert_units(sens, vol.unit); 
% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% conpute an example leadfield
lf = ft_compute_leadfield([0 0 50], sens, vol);

% figure;
% subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf(:,1))
% subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf(:,2))
% subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf(:,3))

