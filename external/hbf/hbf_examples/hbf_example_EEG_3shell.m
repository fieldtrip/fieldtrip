% Short example on the use of Helsinki BEM Framework kernel:
% 3-shell EEG model.
% Please check hbf_example_EMEG_3shell_thorough.m first.
%
% v200928 (c) Matti Stenroos

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
Nmeshes = length(bmeshes);

elecs = hbf_ProjectElectrodesToScalp(electrodes.p_orig,bmeshes);

ci=[1 1/50 1]*.33; 
co=[ci(2:3) 0]; % remember, meshes are ordered in -> out

% make BEM model
D = hbf_BEMOperatorsPhi_LC(bmeshes);
Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,1);
Tphi_elecs = hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,elecs);

% make lead field matrices
LFM_Edir = hbf_LFM_Phi_LC(bmeshes,Tphi_elecs,cortex.p,cortex.nn);
LFM_Exyz = hbf_LFM_Phi_LC(bmeshes,Tphi_elecs,cortex.p);

