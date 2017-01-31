function failed_old_halfspace_bug243

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_old_halfspace_bug243

% this tests for bug 243

% create vol conductor for halfspace (in mm)
geom.pnt = [[0 2 0];[0 0 0];[1 2 0]];
vol_hs = ft_headmodel_halfspace(geom,[-1 -1 -1],'conductivity',1);

% FIXME: this solution takes into account only dipole projected on surface
%   choose adequate electrodes set

  % create vol conductor for sphere (radius = 100 m)
  r = 1e5; % in mm
  pnt = r*icosahedron162;
  pnt(:,3) = pnt(:,3) - r;
  % figure,ft_plot_mesh(pnt),axis on
  vol_sph = ft_headmodel_singlesphere(pnt,'conductivity',1);


% FIXME: test infinite vs halfspace when both electrodes and dipoles are
% very far away from plane

% FIXME: test halfpsace versus BEM when the surface is very big

% create infinite volume conductor
vol = [];
vol.type = 'infinite';

% create electrodes positions
elec.pnt = [[0 0 -0.2];[0 0 -0.4];[0 0 -0.6];[0 0 -0.8];[0 0 -1];[0 0 -1.2];[0 0 -1.4];[0 0 -1.6]];
for i=1:8
  elec.label{i} = sprintf('ch%0.2d',i);
end

% apply the necessary sensors/volume transformations
[vol_hs1, sens1]  = ft_prepare_vol_sens(vol_hs, elec);
% [vol_sph1, sens2] = ft_prepare_vol_sens(vol_sph, elec); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 first try, calculate lf for certain dipoles positions
pos = [[0 -1 -1];[0 -1 -2];[-1 -2.5 -1]];Ndip = size(pos,1);
ori = [[0 0 1];[0 0 1];[0 0 1]]; ori = ori'; mori = ori(:);
[lf]  = ft_compute_leadfield(pos, sens1, vol);  
lf  = lf *mori;
[lf1] = ft_compute_leadfield(pos, sens1, vol_hs1); 
lf1 = lf1 *mori;

figure,plot(lf)
figure,plot(lf1)

% 2 near to surface
pos = [[0 -1 -0.1];[0 -1 -0.2];[-1 -2.5 -0.5]];Ndip = size(pos,1);
[lf]  = ft_compute_leadfield(pos, sens1, vol); 
[lf1] = ft_compute_leadfield(pos, sens1, vol_hs1);
figure,plot(lf)
figure,plot(lf1)

% 3 far away from surface
pos = [[0 -1 -10];[0 -1 -2];[-1 -2.5 -5]];Ndip = size(pos,1);
[lf]  = ft_compute_leadfield(pos, sens1, vol); lf = lf *mori;
[lf1] = ft_compute_leadfield(pos, sens1, vol_hs1);lf1 = lf1 *mori;
figure,plot(lf) 
figure,plot(lf1)

% 4 far away and electrode-dipoles are near
elec.pnt = 1e6*[[0 0.1 -0.2];[0 0 -0.3];[0 0.2 -0.2];[0 0.3 -0.3];[0 0.7 -.2];[0 1 -.3];[0 1.1 -.2];[0 1.4 -.3]];
pos = 1e7*[[0 -1 -1];[0 -1 -2];[0 -2.5 -2]];Ndip = size(pos,1);
[lf]  = ft_compute_leadfield(pos, elec, vol); %lf = lf *mori;
[lf1] = ft_compute_leadfield(pos, elec, vol_hs1); %lf1 = lf1 *mori;
figure,plot(lf) 
figure,plot(lf1)

