function test_headmodel_openmeeg

% MEM 1500mb
% WALLTIME 00:20:00

% TEST ft_prepare_vol_sens ft_compute_leadfield ft_headmodel_openmeeg

% load openmeeg paths
system('module load openmeeg');

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

vol = ft_headmodel_openmeeg(geom, 'conductivity', [1 1/20 1]);

% create a set of electrodes
sel = find(pnt(:,3)>0);
sens.elecpos = pnt(sel,:) * 100;
sens.chanpos = sens.elecpos;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% conpute an example leadfield
lf = ft_compute_leadfield([0 0 50], sens, vol);

if any(mean(lf)./mean(abs(lf)) > 100*eps)
  error('the leadfield is not average reference');
end

figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf(:,3))

vol2 = [];
vol2.o = [0 0 0];
vol2.r = [80 90 100];
vol2.c =[1 1/20 1];
% project the electrodes on the volume conduction model
[vol2, sens2] = ft_prepare_vol_sens(vol2, sens);

% conpute an example leadfield
lf2 = ft_compute_leadfield([0 0 50], sens2, vol2);

if any(mean(lf2)./mean(abs(lf2)) > 100*eps)
  error('the leadfield is not average reference');
end

figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf2(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf2(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf2(:,3))

figure;
subplot(2,2,1); plot(lf2(:,1), lf(:,1), '.')
subplot(2,2,2); plot(lf2(:,2), lf(:,2), '.')
subplot(2,2,3); plot(lf2(:,3), lf(:,3), '.')

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% repeat the test, comparing to a single sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create the BEM geometries
geom = [];
geom.bnd(1).pnt = pnt * 100;
geom.bnd(1).tri = tri;

vol = ft_headmodel_openmeeg(geom, 'conductivity', [1/80]);

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% conpute an example leadfield
lf = ft_compute_leadfield([0 0 50], sens, vol);

if any(mean(lf)./mean(abs(lf)) > 100*eps)
  error('the leadfield is not average reference');
end

figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf(:,3))

vol2 = [];
vol2.o = [0 0 0];
vol2.r = [100];
vol2.c =[1/80];
% project the electrodes on the volume conduction model
[vol2, sens2] = ft_prepare_vol_sens(vol2, sens);

% conpute an example leadfield
lf2 = ft_compute_leadfield([0 0 50], sens2, vol2);

if any(mean(lf2)./mean(abs(lf2)) > 100*eps)
  error('the leadfield is not average reference');
end

figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf2(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf2(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf2(:,3))

figure;
subplot(2,2,1); plot(lf2(:,1), lf(:,1), '.')
subplot(2,2,2); plot(lf2(:,2), lf(:,2), '.')
subplot(2,2,3); plot(lf2(:,3), lf(:,3), '.')




