function test_headmodel_simbio

% TEST test_headmodel_simbio
% TEST ft_headmodel_fem_simbio ft_prepare_vol_sens ft_compute_leadfield ft_headmodel_concentricspheres

% this function tests that simbio forward model works, comparing the results with a 3 concentric
% spheres model

ft_hastoolbox('simbio', 1);

% create 3 spherical meshes and the corresponding volumes
[pnt, tri] = icosahedron162;

% radiuses and origins are defined in mm
svol(1).o = [0,0,0];
svol(1).r = 30;
svol(1).bnd.pnt = pnt;
svol(1).bnd.tri = tri;

svol(2).o = [0,0,0];
svol(2).r = 50;
svol(2).bnd.pnt = svol(2).r*pnt;
svol(2).bnd.tri = tri;

svol(3).o = [0,0,0];
svol(3).r = 60;
svol(3).bnd.pnt = svol(3).r*pnt;
svol(3).bnd.tri = tri;

% % generate a volume of 3 concentric spheres (works if number of voxels is odd)
% res = 1; % in mm
% for i=3:-1:1
%   tmp2 = zeros(151,151,151);
%   xgrid = -svol(i).r:res:svol(i).r;
%   ygrid = xgrid;
%   zgrid = xgrid;
%   [X, Y, Z]  = ndgrid(xgrid, ygrid, zgrid);
%   pos = [X(:) Y(:) Z(:)];
%   [inside] = bounding_mesh(pos, svol(i).bnd.pnt, svol(i).bnd.tri);
%   l = length(xgrid)
%   c = 76; sel = (l-1)./2; % in voxel
%   tmp = reshape(inside,[l l l]); 
%   tmp2(c-sel:c+sel,c-sel:c+sel,c-sel:c+sel) = tmp;
%   MR{i} = tmp2;
% end
% bkgrnd = zeros(151,151,151);
% bkgrnd = MR{1}+MR{2}+MR{3};
% % save('~crimic/test/SimBio/spheres','bkgrnd')
fprintf('Loading a volume with a number N = %d of compartments ... \n', numel(svol))
tmp = load('~crimic/fieldtrip-dev/forward/test/spheres.mat','bkgrnd');
bkgrnd = tmp.bkgrnd;

% generate volume's external surface (mm)
[pnt, tri] = icosahedron162;
o = [0,0,0];
r = 60;
bnd.pnt = r*pnt;
bnd.tri = tri;

% create a set of electrodes (mm)
clear sens
sel = find(bnd.pnt(:,3)>0);
sens.chanpos = bnd.pnt(sel,:);
sens.elecpos = bnd.pnt(sel,:);
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end
sens.type = 'eeg';

% Generate the head model for FEM simbio (conductivities, wireframe, etc.)
% load('~crimic/test/SimBio/spheres','bkgrnd')
transform = eye(4);
transform(1:3,4) = [-76 -76 -76]'; % voxels to mm
% ATTENTION: simbio wants the coordinates in voxel units,
% conductivity in S/m, sensors in voxels, dipoles in voxels
% FIXME: actually sens shouldnt be required, generate fake sens def inside
% this function?
vol  = ft_headmodel_fem_simbio(bkgrnd,'tissue',{'sph1','sph2','sph3'}, ...
                                      'tissueval',[1 2 3], ...
                                      'tissuecond',[0.022 0.33 0.33], ...
                                      'sens',sens, ...
                                      'transform',transform,'unit','mm'); 

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% compute the lead fields in the output voxels
lf  = ft_compute_leadfield([0 0 30], sens, vol);

% plot things
figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf(:,3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare with concentric sphere model
clear geom
for i=1:3
  geom(i).pnt = svol(i).bnd.pnt;
end
vol2 = ft_headmodel_concentricspheres(geom, 'conductivity', 1e-3*[0.022 0.33 0.33]); %S/mm 

% sensors
[pnt, tri] = icosahedron162;
o = [0,0,0];
r = 60;
bnd.pnt = r*pnt;
bnd.tri = tri;

clear sens2
sel = find(bnd.pnt(:,3)>0);
sens2.chanpos = bnd.pnt(sel,:);
sens2.elecpos = sens2.chanpos;
for i=1:length(sel)
  sens2.label{i} = sprintf('chan%03d', i);
end
sens2.type = 'eeg';

% reproject electrodes on the spherical geometry
[vol2, sens2] = ft_prepare_vol_sens(vol2, sens2);

% compute an example leadfield
lf2 = ft_compute_leadfield([0 0 30], sens2, vol2);

figure;
subplot(2,2,1); ft_plot_topo3d(sens2.chanpos, lf2(:,1))
subplot(2,2,2); ft_plot_topo3d(sens2.chanpos, lf2(:,2))
subplot(2,2,3); ft_plot_topo3d(sens2.chanpos, lf2(:,3))
