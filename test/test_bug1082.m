function test_bug1082

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_vol_sens ft_prepare_headmodel ft_compute_leadfield ft_plot_topo3d ft_headmodel_openmeeg

fprintf('***  DIAGNOSTICAL INFORMATION ***\n');
fprintf('test script is running on host: %s\n', gethostname());

% generate a unit sphere
[pnt, tri] = icosahedron162;

% create the BEM geometries (in mm)
bnd = [];
bnd.pnt = pnt * 100;
bnd.tri = tri; % normals outwards
% bnd.tri = fliplr(tri); % normals inwards

% create a set of electrodes
sel = find(pnt(:,3)>0);
sens = [];
sens.chanpos = pnt(sel,:) * 100;
sens.elecpos = sens.chanpos;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end

pos = [0 0 50];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculate volume conductors

% FIXME: it is unclear which convention is used by the singlesphere
% solution, in terms of normals
tmpcfg = [];
tmpcfg.conductivity = [1];
tmpcfg.method = 'singlesphere';
vol1 = ft_prepare_headmodel(tmpcfg,bnd);
[vol1, sens1] = ft_prepare_vol_sens(vol1, sens);
lf1 = ft_compute_leadfield(pos, sens1, vol1);

tmpcfg = [];
tmpcfg.conductivity = [1];  
tmpcfg.method = 'openmeeg';
vol2 = ft_prepare_headmodel(tmpcfg,bnd);
[vol2, sens2] = ft_prepare_vol_sens(vol2, sens);
lf2 = ft_compute_leadfield(pos, sens2, vol2);

tmpcfg = [];
tmpcfg.conductivity = [1];  
tmpcfg.method = 'dipoli';
vol3 = ft_prepare_headmodel(tmpcfg,bnd);
[vol3, sens3] = ft_prepare_vol_sens(vol3, sens);
lf3 = ft_compute_leadfield(pos, sens3, vol3);

% % this calls the singlesphere model
% volt = [];
% volt.r = 100;
% volt.c = 1;
% [volt, senst] = ft_prepare_vol_sens(volt, sens);
% lft = ft_compute_leadfield(pos, senst, volt);


figure;
subplot(2,2,1); ft_plot_topo3d(sens1.chanpos, lf1(:,1))
subplot(2,2,2); ft_plot_topo3d(sens1.chanpos, lf1(:,2))
subplot(2,2,3); ft_plot_topo3d(sens1.chanpos, lf1(:,3))
colorbar

figure;
subplot(2,2,1); ft_plot_topo3d(sens2.chanpos, lf2(:,1))
subplot(2,2,2); ft_plot_topo3d(sens2.chanpos, lf2(:,2))
subplot(2,2,3); ft_plot_topo3d(sens2.chanpos, lf2(:,3))  
colorbar

figure;
subplot(2,2,1); ft_plot_topo3d(sens3.chanpos, lf3(:,1))
subplot(2,2,2); ft_plot_topo3d(sens3.chanpos, lf3(:,2))
subplot(2,2,3); ft_plot_topo3d(sens3.chanpos, lf3(:,3))  
colorbar

% figure;
% subplot(2,2,1); ft_plot_topo3d(senst.chanpos, lft(:,1))
% subplot(2,2,2); ft_plot_topo3d(senst.chanpos, lft(:,2))
% subplot(2,2,3); ft_plot_topo3d(senst.chanpos, lft(:,3))  
% colorbar
