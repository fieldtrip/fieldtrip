% create 3 spherical meshes and the corresponding volumes

addpath /home/common/matlab/fieldtrip_private/
[pnt, tri] = icosahedron642;

% radiuses and origins are defined in cm
svol(1).o = [0,0,0];
svol(1).r = 3;
svol(1).bnd.pnt = pnt;
svol(1).bnd.tri = tri;

svol(2).o = [0,0,0];
svol(2).r = 5;
svol(2).bnd.pnt = svol(2).r*pnt;
svol(2).bnd.tri = tri;

svol(3).o = [0,0,0];
svol(3).r = 6;
svol(3).bnd.pnt = svol(3).r*pnt;
svol(3).bnd.tri = tri;

res = 0.1; % in cm
for i=3:-1:1
  tmp2 = zeros(151,151,151);
  xgrid = -svol(i).r:res:svol(i).r;
  ygrid = xgrid;
  zgrid = xgrid;
  [X, Y, Z]  = ndgrid(xgrid, ygrid, zgrid);
  pos = [X(:) Y(:) Z(:)];
  [inside] = bounding_mesh(pos, svol(i).bnd.pnt, svol(i).bnd.tri);
  l = length(xgrid)
  c = 76; sel = (l-1)./2; % in voxel
  tmp = reshape(inside,[l l l]); 
  tmp2(c-sel:c+sel,c-sel:c+sel,c-sel:c+sel) = tmp;
  MR{i} = tmp2;
end
bkgrnd = zeros(151,151,151);
bkgrnd = MR{1}+MR{2}+MR{3};
% save('~crimic/test/SimBio/spheres','bkgrnd')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concentric spheres volume 
% load('~crimic/test/SimBio/spheres','bkgrnd')
transform = [0.1 0 0 -7.5; 0 0.1 0 -7.5; 0 0 0.1 -7.5; 0 0 0 1]; % voxels to cm

% generates the head model for FEM simbio (conductivities, wireframe, etc.)
% generate sphere's external surface 
[pnt, tri] = icosahedron162;
o = [0,0,0];
r = 6;
bnd.pnt = r*pnt;
bnd.tri = tri;

% FIXME: take care that the units here should be specified and be
% consistent with the wireframe (ft_write_headshape + conversion)

% FIXME: add vgrid to the forward/private path and commit it
vol  = ft_headmodel_fem_simbio(bkgrnd,'tissue',{'sph1','sph2','sph3'}, ...
                                   'tissueval',[1 2 3],'tissuecond',[0.022 0.33 0.33], ...
                                   'transform',transform,'unit','cm', ...
                                   'bnd',bnd); 

% create a set of electrodes
clear sens
sel = find(bnd.pnt(:,3)>0);
sens.pnt = bnd.pnt(sel,:);
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end
sens.type = 'eeg';

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% compute the lead fields in the output voxels
lf  = ft_compute_leadfield([0 0 3], sens, vol);
figure;
subplot(2,2,1); ft_plot_topo3d(sens.pnt, lf(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.pnt, lf(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.pnt, lf(:,3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare with concentric sphere model
clear geom
for i=1:3
  geom(i).pnt = svol(i).bnd.pnt;
end
vol2 = ft_headmodel_concentricspheres(geom, 'conductivity', [0.022 0.33 0.33]);

% sensors
sel = find(bnd.pnt(:,3)>0);
sens2.pnt = bnd.pnt(sel,:);
for i=1:length(sel)
  sens2.label{i} = sprintf('chan%03d', i);
end
sens2.type = 'eeg';

% % project sens
% [vol2, sens2] = ft_prepare_vol_sens(vol2, sens2);

% conpute an example leadfield
lf2 = ft_compute_leadfield([0 0 3], sens2, vol2);

figure;
subplot(2,2,1); ft_plot_topo3d(sens2.pnt, lf2(:,1))
subplot(2,2,2); ft_plot_topo3d(sens2.pnt, lf2(:,2))
subplot(2,2,3); ft_plot_topo3d(sens2.pnt, lf2(:,3))
