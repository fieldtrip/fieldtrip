function test_ft_read_and_plot_atlas

% WALLTIME 00:15:00
% MEM 8gb
% DEPENDENCY ft_read_atlas ft_sourceplot

% spm8 might have mexfile issues
ft_hastoolbox('spm12',1);
[ftver, ftpath] = ft_version;

% load MNI pial mesh
pial = load(fullfile(ftpath, '/template/anatomy/surface_pial_both.mat'));

% read and plot AFNI
atlas      = ft_read_atlas(fullfile(ftpath, '/template/atlas/afni/TTatlas+tlrc.HEAD'));
atlas.tissue = double(atlas.brick0(:,:,:,1));
atlas.dim = size(atlas.tissue);
atlas.coordsys = 'mni';

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot AAL
atlas      = ft_read_atlas(fullfile(ftpath, '/template/atlas/aal/ROI_MNI_V4.nii'));

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot Brainweb
brainweb = load(fullfile(ftpath, '/template/atlas/brainweb/brainweb_discrete.mat'));
atlas      = brainweb.atlas; clear brainweb

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot JuBrain
atlas = ft_read_atlas(fullfile(ftpath, '/template/atlas/spm_anatomy/AllAreas_v18_MPM.mat'));

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot VTPM
load(fullfile(ftpath, '/template/atlas/vtpm/vtpm.mat'));
atlas = vtpm; clear vtpm

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot Brainnetome
atlas = ft_read_atlas(fullfile(ftpath, '/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii'));

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% Plot Melbourne Subcortical Atlas
atlas = ft_read_atlas(fullfile(ftpath, '/template/atlas/melb_subcortical/melb_sub.mat'));

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether some other AFNI atlases can be read as well
ft_hastoolbox('afni', 1);

t = tempdir;
untar('https://afni.nimh.nih.gov/pub/dist/atlases/afni_atlases_dist.tgz',t);

% test script for updated BrikInfo.m code

% the datadir contains the unzipped afni_atlases_distr.tgz
datadir = fullfile(t, 'afni_atlases_dist');
f = dir(fullfile(datadir, '*.BRIK.gz'));

err = false(numel(f),1);
for k = 1:numel(f)
  fname = fullfile(f(k).folder, f(k).name);
  [err(k,1), info{k}] = BrikInfo(fname);
end

% check that the reading did not yield errors
assert(all(err==0));

% the fieldtrip code in https://github.com/schoffelen/fieldtrip/tree/afni
% should be able to read in all images without error
ok = true;
try
  for k = 1:numel(f)
    fname = fullfile(f(k).folder, f(k).name);
    x{k} = ft_read_mri(fname);
  end
catch
  ok = false;
end
assert(ok==1);

hastable = false(numel(f),1);
for k = 1:numel(info)
  if isfield(info{k}, 'ATLAS_LABEL_TABLE')
    hastable(k) = true;
  end
end

% the fieldtrip code should be able to read in all atlases without error
ok = true;
%try
  for k = find(hastable)'
    fname = fullfile(f(k).folder, f(k).name);
    y{k} = ft_read_atlas(fname);
  end
%catch
%  ok = false;
%end
%assert(ok==1);

