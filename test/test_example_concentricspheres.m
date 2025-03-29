function test_example_concentricspheres

% WALLTIME 01:00:00
% MEM 10gb
% DEPENDENCY ft_prepare_headmodel ft_prepare_mesh ft_datatype_segmentation
% DATA public

%% read mri
mri = ft_read_mri(dccnpath('/project/3031000.02/external/download/test/ctf/Subject01.mri'));

% check that the anatomical MRI is expressed in millimeter
disp(mri.unit)

%% align the MRI with the desired head coordinate system

ft_determine_coordsys(mri, 'interactive', 'no');

% the specific MRI is already aligned with the CTF coordinate system
% but if yours is not, you should use FT_VOLUMEREALIGN to coregister it
% see https://www.fieldtriptoolbox.org/faq/coordsys/

%% segmentation
cfg          = [];
cfg.output   = {'brain', 'skull', 'scalp'};
segmentedmri = ft_volumesegment(cfg, mri);

% check that the segmentation is also expressed in millimeter
disp(segmentedmri.unit)

%% read electrodes and align them with the MRI 

% determine where FieldTrip and the electrode templates are installed
[ftver, ftpath] = ft_version;

% read some template electrodes
elec = ft_read_sens(fullfile(ftpath,'template/electrode/standard_1020.elc'));

% check the units of the electrode positions, these might be expressed in millimeter or centimeter
disp(elec.unit)

% check which electrode positions are specified, we expect these to include
% anatomical landmarks or fiducials
disp(elec.label)

% we will use the fiducials that were identified in the anatomical MRI for coregistration
% in this case the fiducials are stored in the mri.fid field
nas = mri.fid.pos(1,:);
lpa = mri.fid.pos(2,:);
rpa = mri.fid.pos(3,:);

% since the fiducials from the MRI are still expressed in millimeter, the electrodes should be as well
elec = ft_convert_units(elec, 'mm');

% the mri fiducials are named 'nas/lpa/rpa', but the electrodes include 'Nz/LPA/RPA'
fiducials.pos     = [nas; lpa; rpa];
fiducials.label   = {'Nz','LPA','RPA'};
fiducials.unit    = 'mm';

% match the fiducials from the electrodes with those of the MRI
cfg          = [];
cfg.method   = 'fiducial';
cfg.target   = fiducials;
cfg.elec     = elec;
cfg.fiducial = {'Nz', 'LPA', 'RPA'}; % these labels should match with both the fiducials and elec structure
elec_aligned = ft_electroderealign(cfg);

% the ear points in the MRI are placed at the ear canals, but at the pre-auricular points in the electrodes
% an interactive realignment is recommended to check that it all fits
for i=1:size(elec_aligned.chanpos,1)
  % add 12 mm to the x-axis, needed due to different definitions of the ear points
  elec_aligned.chanpos(i,1) = elec_aligned.chanpos(i,1) + 12;
  elec_aligned.elecpos(i,1) = elec_aligned.elecpos(i,1) + 12;
end

%% make the volume conductor

% convert the segmented MRI from millimeter to meter
segmentedmri = ft_convert_units(segmentedmri, 'm');
disp(segmentedmri.unit)

% convert the electrodes from millimeter to meter
elec_aligned = ft_convert_units(elec_aligned, 'm');
disp(elec_aligned.unit)

cfg              = [];
cfg.method       = 'concentricspheres';
cfg.tissue       = {'brain', 'skull', 'scalp'};
cfg.conductivity = [0.33 0.0125 0.33]; % this is Siemens per meter
headmodel        = ft_prepare_headmodel(cfg, segmentedmri);

%% make a regular 3D grid as the sourcemodel

cfg           = [];
cfg.method    = 'basedongrid';
cfg.xgrid     = -10:1:10;
cfg.ygrid     = -10:1:10; % from -10 to 10 in steps of 1
cfg.zgrid     =  -2:1:10;
cfg.unit      = 'cm';
cfg.tight     = 'yes';
cfg.headmodel = headmodel; % this determines which dipoles are inside/outside
sourcemodel = ft_prepare_sourcemodel(cfg);

% convert the sourcemodel from centimeter to meter
sourcemodel = ft_convert_units(sourcemodel, 'm');

% plot the headmodel together with the source positions

% note that there are three spheres in the figure below 
% and that the electrodes are ideally to be on the surface of the outermost (scalp) sphere

figure
ft_plot_headmodel(headmodel, 'facecolor', 'skin')
ft_plot_sens(elec_aligned, 'elecsize', 0.010, 'elecshape', 'disc', 'facecolor', 'k', 'label', 'label')
ft_plot_mesh(sourcemodel.pos)
ft_plot_axes([], 'unit', 'm', 'coordsys', 'ctf')
ft_headlight
alpha 0.5

%% compute the leadfields

% take a subset of 29 channels
chansel = ft_channelselection('eeg1020', elec_aligned);

% remove some that we are not interested in
chansel = setdiff(chansel, 'T3'); % this is an alternative name for C7, so it is double
chansel = setdiff(chansel, 'T4'); % this is an alternative name for C8, so it is double
chansel = setdiff(chansel, 'T5'); % this is an alternative name for P7, so it is double
chansel = setdiff(chansel, 'T6'); % this is an alternative name for P8, so it is double
chansel = setdiff(chansel, 'A1'); % left earlobe
chansel = setdiff(chansel, 'A2'); % right earlobe

cfg             = [];
cfg.headmodel   = headmodel;    % expressed in meter and Siemens/meter
cfg.elec        = elec_aligned; % expressed in meter
cfg.sourcemodel = sourcemodel;  % expressed in meter
cfg.channel     = chansel;
leadfield       = ft_prepare_leadfield(cfg); % expressed in Volt per Ampere*meter

%% plot the topography for a given source

dippos = [0 0 50]/1000;     % in meter
dipori = [1 0 0]';          % along the x-direction 
dipmom = 100*1e-9 * dipori; % 100 nAm

% find the grid position that is the closest to the desired dipole position
[~, dipindx] = min(pdist2(leadfield.pos, dippos));

% show the position that was found
disp(leadfield.pos(dipindx,:));

% compute the potential distribution
potential = leadfield.leadfield{dipindx} * dipmom;

% make an electrode structure with the same subset
elecindx = match_str(elec_aligned.label, chansel);

elec1020 = [];
elec1020.label = elec_aligned.label(elecindx);
elec1020.elecpos = elec_aligned.elecpos(elecindx,:);
elec1020.chanpos = elec_aligned.chanpos(elecindx,:);

figure
% ft_plot_headmodel(headmodel, 'facecolor', 'skin')
ft_plot_sens(elec1020, 'elecsize', 0.010, 'elecshape', 'disc', 'facecolor', 'k', 'label', 'label')
ft_plot_dipole(dippos, dipori, 'unit', 'm')
ft_plot_topo3d(elec1020.elecpos, potential)
ft_plot_axes([], 'unit', 'm', 'coordsys', 'ctf')
ft_headlight
alpha 0.8
colorbar
