function test_pullXXX

% MEM gb
% WALLTIME
% DEPENDENCY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function creates a tetrahedral two-sphere mesh with cylindrical electrodes lying on the boundary between the two spheres (i.e., intracranial electrodes)
% which can be used to perform (e.g., ECoG) forward simulations, particularly in FEMfuns (Finite Element Method for useful neuroscience simulations, https://github.com/meronvermaas/FEMfuns).
% It has been created in the IntoTheBrain project (@mcpiastra @ThomOostendorp) and functions will be extended to allow the full use of FEMfuns within FieldTrip.
% The structure of this script is more or less:
% 1. Create two triangulated spheres representing the brain and skull and electrode coordinates on top of the brain surface
% 2. Create cylinders at each electrode coordinate.
%   a. make the cylinder surface crossing through the spherical surface it is combined with
%   b. refine the spherical surface where the cylinder crosses it
%   c. cut off the part of the cylindrical surface sticking into the sphere
% 3. Combine all the surfaces into one surface with unique faces and nodes
% 4. Make a tetrahedral mesh

% this requires the external iso2mesh toolbox
ft_hastoolbox('iso2mesh', 1);

%create spheres representing the brain and skull
c0 = [0,0,0];
%'brain' sphere
rad_brain = 7;
[brainsurf.node,brainsurf.face]=meshasphere(c0,rad_brain,0.3,3);
%'skull' sphere
rad_skull = 10;
[skullsurf.node,skullsurf.face]=meshasphere(c0,rad_skull,0.3,3);
%get labeling points of the brain and skull surface
bpnt = point_in_surf(brainsurf);
spnt = point_in_surf(skullsurf);

% define some electrode coordinates
N = 3;
theta=linspace(0,2*pi,N);
phi=linspace(0,pi,N);
[theta,phi]=meshgrid(theta,phi);

x = rad_brain*sin(phi).*cos(theta);
y = rad_brain*sin(phi).*sin(theta);
z = rad_brain*cos(phi);
elec_coords = [x(:),y(:),z(:)];
elec_coords = unique(round(elec_coords,2),'rows');

% combine the inner skull surface (i.e., brain sphere) with electrode surfaces and get labeling points of the electrodes
dp_elec = 0.5; %height  of the electrode cylinder
rad_elec = 0.2; %radius of the electrode cylinder
[skullelsurf,elecmarkers] = add_electrodes(brainsurf, elec_coords, rad_elec, dp_elec); %TODO: write script to check consistency of normal directions (in vs out)
% plotmesh(skullelsurf.node, skullelsurf.face) %inspect the added electrodes

%combine with all other surfaces; in this case only with outer skull sphere
%slow but correct; TODO compile the add_surf function or vectorize to increase speed
merged_surfs = add_surf(skullsurf,skullelsurf);

%create volumetric tetrahedral mesh
[tet_node,tet_elem,tet_face] = s2m(merged_surfs.node,merged_surfs.face, 1, 1, 'tetgen', [bpnt;spnt;elecmarkers]);
%[tet_node,tet_elem,tet_face] = s2m(merged_surfs.node,merged_surfs.face, 1, 1, 'tetgen1.5'); %automatic labeling, order of labels unclear

%save the tetrahedral mesh for inspection
%savemsh(tet_node,tet_elem,'spheres_elecs.msh')
