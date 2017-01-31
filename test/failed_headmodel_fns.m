function failed_headmodel_fns

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_ft_compute_leadfield
% TEST ft_compute_leadfield ft_headmodel_fns ft_prepare_vol_sens ft_compute_leadfield

% this function tests that FNS forward model works, comparing the results with a 3 concentric
% spheres model

ft_hastoolbox('fns', 1);

% create 3 spherical meshes and the corresponding volumes
[pnt, tri] = icosahedron162;

% radiuses and origins are defined in mm
svol(1).o = [0,0,0];
svol(1).r = 10;
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

% load the correspondent volumetric matrix
fprintf('Loading a volume with a number N = %d of compartments ... \n', numel(svol));
tmp = load(dccnpath('/home/common/matlab/fieldtrip/data/test/spheres.mat'));
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

% Generate the head model for FNS (conductivities, labels, etc.)
% load('~crimic/test/SimBio/spheres','bkgrnd')
transform = eye(4);
transform(1:3,4) = [-76 -76 -76]'; % voxels to mm

% FIXME: write a function that reads in the segmentation fields and
% transforms them into a integer compartments image

% ATTENTION: fns wants the coordinates in voxel units,
% conductivity in S/m, sensors in voxels, dipoles in voxels
vol  = ft_headmodel_fns(bkgrnd,'tissue',{'sph1','sph2','sph3'}, ...
                                      'tissueval',[1 2 3], ...
                                      'tissuecond',[0.022 0.33 0.33], ...
                                      'sens',sens, ...
                                      'transform',transform,'unit','mm'); 

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);
                    
% compute the lead fields in the output voxels
lf  = ft_compute_leadfield([0 0 30], sens, vol); % actually here sens not necessary anymore

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
