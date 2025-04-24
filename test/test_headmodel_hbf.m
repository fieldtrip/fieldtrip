function test_headmodel_hbf

% MEM 1gb
% WALLTIME 00:20:00
% DEPENDENCY ft_prepare_vol_sens ft_compute_leadfield ft_headmodel_openmeegf

% DATA no

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG, comparing to 3 sphere model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate a unit sphere
[pos, tri] = mesh_sphere(162);

% create the BEM geometries
bnd = [];
bnd(1).pos = pos * 0.1;
bnd(1).tri = tri;
bnd(1).unit = 'm';
bnd(2).pos = pos * 0.09;
bnd(2).tri = tri;
bnd(2).unit = 'm';
bnd(3).pos = pos * 0.08;
bnd(3).tri = tri;
bnd(3).unit = 'm';

vol = ft_headmodel_hbf(bnd, 'conductivity', [1 1/20 1; 1/20 1 0]);

% create a set of electrodes
sel = find(pos(:,3)>0);
sens.elecpos = pos(sel,:) * 0.1;
sens.chanpos = sens.elecpos;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% conpute an example leadfield
lf = ft_compute_leadfield([0 0 0.05], sens, vol);

if any(mean(lf)./mean(abs(lf)) > 100*eps)
  error('the leadfield is not average reference');
end

figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf(:,3))

vol2 = [];
vol2.o = [0 0 0];
vol2.r = [0.08 0.09 0.1];
vol2.c =[1 1/20 1];
vol2.unit = 'm';
% project the electrodes on the volume conduction model
[vol2, sens2] = ft_prepare_vol_sens(vol2, sens);


% conpute an example leadfield
lf2 = ft_compute_leadfield([0 0 0.05], sens2, vol2);

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
% MEG, comparing to a single sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create the BEM geometries
bnd = [];
bnd(1).pos = pos * 0.1;
bnd(1).tri = tri;
bnd(1).unit = 'm';

vol = ft_headmodel_hbf(bnd, 'conductivity', [1/80; 0]);

% create some coils
sel = find(pos(:,3)>0);
sens = [];
sens.coilpos = pos(sel,:) * 0.12;
sens.coilori = pos(sel,:);
sens.chanpos = sens.coilpos;
sens.chanori = sens.coilori;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% conpute an example leadfield
lf = ft_compute_leadfield([0 0 0.05], sens, vol, 'reducerank', 2);

if any(lf(:,3) > eps)
  error('silent compoenent is not silent!');
end

figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf(:,3))

vol2 = [];
vol2.o = 0;
vol2.r = 0.1;
% project the electrodes on the volume conduction model
[vol2, sens2] = ft_prepare_vol_sens(vol2, sens);

% conpute an example leadfield
lf2 = ft_compute_leadfield([0 0 0.05], sens2, vol2, 'reducerank', 2);
vecnorm(lf2)

figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf2(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf2(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf2(:,3))

figure;
subplot(2,2,1); plot(lf2(:,1), lf(:,1), '.')
subplot(2,2,2); plot(lf2(:,2), lf(:,2), '.')
subplot(2,2,3); plot(lf2(:,3), lf(:,3), '.')
