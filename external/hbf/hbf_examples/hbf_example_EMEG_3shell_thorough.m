% Example on the use of Helsinki BEM Framework kernel:
% 3-shell MEG+EEG model.
% 
% Please start by downloading the sample data from
% github.com/MattiStenroos/hbf_Sampledata
%
% v200924 (c) Matti Stenroos

% The example uses the head geometry of Subject 1 of the data set described in
% Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal
% human neuroimaging dataset. Sci. Data 2:150001 doi: 10.1038/sdata.2015.1,
% available through 
% ftp://ftp.mrc-cbu.cam.ac.uk/personal/rik.henson/wakemandg_hensonrn/
%
% The meshes were constructed with SPM, using scripts provided with the
% data. If you use the head geometry in a publication,
% please cite the Wakeman--Henson publication and include the above
% information.
%
% Sensor geometries of a 306-channel Neuromag MEG system (as described in
% MNE-C software, http://mne.tools/) and 256-channel ABC EEG-layout (as
% described by Biosemi, http://www.biosemi.com/download/Cap_coords_all.xls,
% were added manually by the author. The sensors and their coregistrations
% do not correspond to the sensors used in the of the Wakeman--Henson
% dataset.
%
% Please download the head geometry data from
% github.com/MattiStenroos/hbf_sampledata

clear
% Add BEM framework paths. Put here your hbf root directory
hbfroot = '~/gitwork/hbf_lc_open/';
thisdir = cd;
cd(hbfroot);
hbf_SetPaths();
cd(thisdir);

% put here your path to the example data
sampledatapath = '~/gitwork/hbf_sampledata/hbf_samplehead_3shell_wh';
load(sampledatapath);
% -> bmeshes, cortex, electrodes, coils

% 1. Boundary meshes
%
% Boundary meshes are described as a Nx1,cell array  where each cell is a
% struct that contains fields "p" for points (vertices), and "e" for
% element description (faces). The element description is 1-based.
% bmeshes{I} =
%
%   p: [Number of vertices x 3]
%   e: [Number of triangles x 3]
%
%
% First time with new meshing: tests...
% Problems with the BEM are typically due to unsuitable or ill-specified
% meshes. Even if meshes are otherwise good, the orientation of the
% triangles is often wrong; this is the most common reason for computations
% going wrong. This BEM framework assumes CCW orientation. To check and
% flip the orientation, you can just type

Nmeshes = length(bmeshes);
success = zeros(Nmeshes,1);
for M = 1:Nmeshes
    [bmeshes{M},success(M)] = hbf_CorrectTriangleOrientation(bmeshes{M});
end

% You can do further checks using, e.g., this simple tool. If orientation
% correction failed, this may give you at least some tips, where to look
% for the error.

status = cell(Nmeshes,1);
for M = 1:Nmeshes
    status{M} = hbf_CheckMesh(bmeshes{M});
end

% The convention of this BEM framework is that the innermost mesh has number
% 1 and outermost mesh M; this is mandatory. If the model is nested (like a
% 3-shell model), the order can be checked/set with
bmeshes = hbf_SortNestedMeshes(bmeshes);



% 2. MEG sensors
%
% Set MEG coils. There are some options for describing coils:
% coils = 
%   QP: field computations points [Number of field computation points x 3]
%   QN: sensor orientation        [Number of field computation points x 3];
%       each must be unit length
%   QPinds: start and end indices of the QPs of each sensor [number of sensors x 2]
%   QW: sensor integration weights [number of field computation points x 1]
%   QtoC: integral weights and conversion QP -> sensors
%      [Number of sensors x Number of field computation points]
%   p:  sensor positions [Number of sensors x 3]
%   n:  sensor normals [Number of sensors x 3]; each must be unit-length.
%
% You must give (QP, QN and QtoC) or (QP, QN, QP and QPinds) or (p and n).
% If overlapping or conflicting information is given, the order of
% preference is the same as in this list.
% 
% This sample coil set represents a 306-channel Neuromag system and uses
% 4-point sensor integrals as in, for example, Neuromag and MNE software. The
% coregistration is arbitrary, not the one in the data set. There is some
% further information in fields 'names','mesh3d', and 'description'.

% 3. EEG electrodes
%
% Set of EEG electrodes.

% The electrodes are assumed to lie on scalp; the given electrode points
% need to be roughly coregistered to scalp mesh. If you wish to place
% electrodes elsewhere, contact me. In this case, we have 256 electrodes.
% For convenience, the imported struct 'electrodes' contains, in addition
% to locations, also electrode names and a 3D mesh (unprojected positions).

% Project electrodes to scalp and build interpolation operator from model
% nodes to projected electrodes. This information will be collected in 
% 'elecs' struct.
elecs = hbf_ProjectElectrodesToScalp(electrodes.p_orig,bmeshes);
% let's include also these, although these are not needed for computations 
elecs.names = electrodes.names;
elecs.mesh3d = electrodes.mesh3d;
elecs.mesh3d.p = elecs.p_proj;

% 4. Source space
% We need to provide dipole locations. If we want to use oriented sources,
% we need to provide also the orientations.
%
% Here we have a cortical mesh with also normal vectors for mesh vertices>
% cortex = 
% 
%   struct with fields:
% 
%      p: [20484×3 single]
%      e: [40960×3 int32]
%     nn: [20484×3 double]
%
% so we can directly use cortex.p as source positions and cortex. nn as
% source locations.

% 5. Conductivities
% Set conductivities for the head model. In this example, we use a (nested)
% 3-layer model, so we need to give conductivities for the brain, skull,
% and scalp compartments. Remember to give the conductivities in the
% correct order!

ci=[1 1/50 1]*.33; % conductivity inside each surface --- remember the order!
co=[ci(2:3) 0];    % conductivity outside each surface

% ------------------------------------------------------------------------
% Now we should be ready for computation, unless something went wrong with
% the checks...

starttime = clock;

% 1. BEM double-layer operators for potentials: do this once per meshing...
D = hbf_BEMOperatorsPhi_LC(bmeshes);

% 2. BEM operators for magnetic field due to volume currents: 
% do this, if the boundary meshing or MEG sensors change

DB = hbf_BEMOperatorsB_Linear(bmeshes,coils);

%3. Full transfer matrix
% Build BEM transfer matrix for potential using the isolated source
% approach, isolation set to surface 1.

Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,1);

%4. Transfer matrices for volume component of the magnetic field and the
%potential

TBvol = hbf_TM_Bvol_Linear(DB,Tphi_full,ci,co);
% Extract/interpolate EEG transfer matrix for electrodes only
Tphi_elecs = hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,elecs);


%% 5. Forward solutions = leadfield matrices
% Compute directed LFM --- source directions must have unit norm!
% You can use this function as well for computing the magnetic field due to
% any set of dipoles; in that case, just use the dipole moment vectors
% instead of unit-length source direction vectors.

[LFM_Edir, LFM_Mdir] = hbf_LFM_PhiAndB_LC(bmeshes,coils,Tphi_elecs,TBvol,cortex.p,cortex.nn);

% LFM for orthogonal unit dipoles (orientations x,y,z in world coordinates).

[LFM_Exyz, LFM_Mxyz] = hbf_LFM_PhiAndB_LC(bmeshes,coils,Tphi_elecs,TBvol,cortex.p);
fprintf('Total time was %ds.\n',round(etime(clock,starttime)));

% % MEG only...
% LFM_Mdir = hbf_LFM_B_LC(bmeshes,coils,TBvol,cortex.p,cortex.nn);
% LFM_Mxyz = hbf_LFM_B_LC(bmeshes,coils,TBvol,cortex.p);

% % EEG only...
% LFM_Edir = hbf_LFM_Phi_LC(bmeshes,Tphi_elecs,cortex.p,cortex.nn);
% LFM_Exyz = hbf_LFM_Phi_LC(bmeshes,Tphi_elecs,cortex.p);




