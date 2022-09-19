function test_example_fem

% WALLTIME 08:00:00
% MEM 12gb
% DEPENDENCY ft_prepare_headmodel ft_prepare_mesh ft_datatype_segmentation

[ftver, ftpath] = ft_version;

%% read mri
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.mri'));

%% segmentation
cfg          = [];
cfg.output   = {'gray', 'white', 'csf', 'skull', 'scalp'};
segmentedmri = ft_volumesegment(cfg, mri);

%% mesh
cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
mesh       = ft_prepare_mesh(cfg,segmentedmri);

%% volume conductor
cfg              = [];
cfg.method       = 'simbio';
cfg.conductivity = [0.33 0.14 1.79 0.01 0.43];
headmodel        = ft_prepare_headmodel(cfg, mesh);

%% electrode alignment
elec = ft_read_sens(fullfile(ftpath,'template/electrode/standard_1020.elc'));

nas = mri.hdr.fiducial.head.nas;
lpa = mri.hdr.fiducial.head.lpa;
rpa = mri.hdr.fiducial.head.rpa;

%
fiducials.pos     = [nas; lpa; rpa];
fiducials.label   = {'Nz','LPA','RPA'};
fiducials.unit    = 'mm';

cfg          = [];
cfg.method   = 'fiducial';
cfg.target   = fiducials;
cfg.elec     = elec;
cfg.fiducial = {'Nz', 'LPA', 'RPA'};
elec_align   = ft_electroderealign(cfg);

% add 12 mm to x-axis
n=size(elec_align.chanpos,1);
for i=1:n
 elec_align.chanpos(i,1)=elec_align.chanpos(i,1)+12;
 elec_align.elecpos(i,1)=elec_align.elecpos(i,1)+12;
end

%% make the sourcemodel/grid

% At the moment the sourcemodel is defined prior
% to the leadfield because ft_prepare_sourcemodel does not automatically create
% a sourcemodel based on a hexahedral headmodel.

cfg                 = [];
cfg.mri             = mri;
cfg.resolution      = 3;
sourcemodel         = ft_prepare_sourcemodel(cfg);
sourcemodel         = ft_convert_units(sourcemodel, headmodel.unit);
sourcemodel.inside(500:end) = false; % do only a few dipoles

cfg            = [];
cfg.headmodel  = headmodel;
cfg.elec       = elec_align;
cfg.sourcemodel = sourcemodel;
cfg.channel    = elec_align.label(4:6); % do only a few channels, to save time
lf             = ft_prepare_leadfield(cfg);

m = cellfun(@mean,lf.leadfield(lf.inside),'uniformoutput',false)';
m = cat(1, m{:});
assert(all(m(:)<eps));
