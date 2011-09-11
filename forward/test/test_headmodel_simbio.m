% generate a unit sphere
[pnt, tri] = icosahedron162;

% generate sphere external surface (radiuses: 3,5,6 cm)
o = [0,0,0];
r = 6;
bnd.pnt = r*pnt;
bnd.tri = tri;

% concentric spheres volume (~/test/SimBio/test_algo.m)
load('~crimic/test/SimBio/spheres','bkgrnd')
transform = [0.1 0 0 -7.5; 0 0.1 0 -7.5; 0 0 0.1 -7.5; 0 0 0 1]; % voxels to cm

% generates the head model for FNS (conductivities and segmentation)
% FIXME: do i need a bnd to represent the surface of the head here?

% FIXME: take care that the units here should be specified and be
% consistent with the wireframe (ft_write_headshape + conversion)

% FIXME: add vgrid to the forward/private path and commit it
vol  = ft_headmodel_fem_simbio(bkgrnd,'tissue',{'sph1','sph2','sph3'}, ...
                                   'tissueval',[1 2 3],'tissuecond',[0.022 0.33 0.33], ...
                                   'transform',transform,'unit','cm', ...
                                   'bnd',bnd); 

% create a set of electrodes
sel = find(pnt(:,3)>0);
sens.pnt = pnt(sel,:) * 6;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end
sens.type = 'eeg';

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);

% compute the lead fields in the output voxels
lf  = ft_compute_leadfield([0 0 3], sens, vol);
