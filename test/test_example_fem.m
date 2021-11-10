function test_example_fem

% WALLTIME 08:00:00
% MEM 12gb
% DEPENDENCY ft_prepare_headmodel ft_prepare_mesh ft_datatype_segmentation

[ftver, ftpath] = ft_version;

%% read mri
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));

%% segmentation
cfg          = [];
cfg.output   = {'gray', 'white', 'csf', 'skull', 'scalp'};
segmentedmri = ft_headmodelumesegment(cfg, mri);

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
sourcemodel         = ft_prepare_sourcemodel(cfg);
sourcemodel         = ft_convert_units(sourcemodel, headmodel.unit);

cfg            = [];
cfg.headmodel  = headmodel;
cfg.elec       = elec_align;
cfg.sourcemodel = sourcemodel;
lf             = ft_prepare_leadfield(cfg);

% plot the leadfield for a few representative locations: points around z-axis with increasing z values

plotpos = [];
n=size(lf.pos,1);
p=1;
for i = 1:n
  if lf.pos(i,1)==-0.1 && lf.pos(i,2)==-0.2
      plotpos(p)=i;
      p=p+1;
  end
end

figure;
for i=1:20

  subplot(4,5,i);
  ft_plot_topo3d(lf.cfg.elec.chanpos,lf.leadfield{plotpos(i)}(:,3));
  % view([0 0]);

end
