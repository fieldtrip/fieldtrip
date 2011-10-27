% generate a unit sphere
[pnt, tri] = icosahedron162;

% create the BEM geometries (in mm)
bnd = [];
bnd.pnt = pnt * 100;
bnd.tri = tri;

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
tmpcfg = [];
tmpcfg.conductivity = [1];
tmpcfg.method = 'singlesphere';
vol1 = ft_prepare_headmodel(tmpcfg,bnd);
[vol1, sens1] = ft_prepare_vol_sens(vol1, sens);
lf1 = ft_compute_leadfield(pos, sens1, vol1);

tmpcfg = [];
tmpcfg.conductivity = [1];  
tmpcfg.method = 'bem_openmeeg';
vol2 = ft_prepare_headmodel(tmpcfg,bnd);
[vol2, sens2] = ft_prepare_vol_sens(vol2, sens);
lf2 = ft_compute_leadfield(pos, sens, vol2);

% this calls the singlesphere model
volt = [];
volt.r = 100;
volt.c = 1;
[volt, senst] = ft_prepare_vol_sens(volt, sens);
lft = ft_compute_leadfield(pos, senst, volt);


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
subplot(2,2,1); ft_plot_topo3d(senst.chanpos, lft(:,1))
subplot(2,2,2); ft_plot_topo3d(senst.chanpos, lft(:,2))
subplot(2,2,3); ft_plot_topo3d(senst.chanpos, lft(:,3))  
colorbar
