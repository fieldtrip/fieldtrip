function failed_bug2396

% WALLTIME 00:10:00
% MEM 2000mb

% TEST test_bug2396
% TEST ft_read_headshape

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2396'));

%% define the files

mrifile = 'nobias_KR.nii';
labelfile = 'segmentation/brainsuite/nobias_KR.label.nii';

dfsfile = {
  'segmentation/brainsuite/nobias_KR.brain.dfs'
  'segmentation/brainsuite/nobias_KR.inner_skull.dfs'
  'segmentation/brainsuite/nobias_KR.outer_skull.dfs'
  'segmentation/brainsuite/nobias_KR.scalp.dfs'
  };

meshfile = {
  'segmentation/brainvisa/nobias_KR_BEM_brain.mesh'
  'segmentation/brainvisa/nobias_KR_BEM_head.mesh'
  'segmentation/brainvisa/nobias_KR_BEM_skull.mesh'
  'segmentation/brainvisa/nobias_KR_head.mesh'
  };

% the following contain metadata that goes with the mesh files
% segmentation/brainvisa/nobias_KR_BEM_brain.mesh.minf
% segmentation/brainvisa/nobias_KR_BEM_head.mesh.minf
% segmentation/brainvisa/nobias_KR_BEM_skull.mesh.minf
% segmentation/brainvisa/nobias_KR_head.mesh.minf

surffile = {
  'segmentation/freesurfer/lh.BEM_brain_surface.surf'
  'segmentation/freesurfer/lh.BEM_inner_skull_surface.surf'
  'segmentation/freesurfer/lh.BEM_outer_skin_surface.surf'
  'segmentation/freesurfer/lh.BEM_outer_skull_surface.surf'
  };

%% read the content from the files

mri = ft_read_mri(mrifile);

% it would be more consistent to read it with ft_read_atlas, as if it were a ft_datatype_segmentation
seg = ft_read_mri(labelfile);

% dress up the segmentation
seg.anatomylabel = {};
seg.anatomylabel{19} = 'brain';
seg.anatomylabel{18} = 'csf';
seg.anatomylabel{17} = 'skull';
seg.anatomylabel{16} = 'scalp';
seg.coordsys = 'ctf'; % it seems to be in CTF coordinates, i.e. positive y axis through nasion

cfg = [];
cfg.atlas = seg; % use itself as the atlas to look up the number-name mappings
ft_sourceplot(cfg, seg);

%% read the content from the files

dfs = {};
figure
for i=1:4
  dfs{i} = ft_read_headshape(dfsfile{i});
  ft_plot_mesh(dfs{i}, 'edgecolor', 'none', 'facecolor', [1 1 1]/2, 'facealpha', 0.3);
end

mesh = {};
figure
for i=1:4
  mesh{i} = ft_read_headshape(meshfile{i});
  ft_plot_mesh(mesh{i}, 'edgecolor', 'none', 'facecolor', [1 1 1]/2, 'facealpha', 0.3);
  % note that mesh number 4 is ugly and has internal structures
end


surf = {};
figure
for i=1:4
  surf{i} = ft_read_headshape(surffile{i});
  ft_plot_mesh(surf{i}, 'edgecolor', 'none', 'facecolor', [1 1 1]/2, 'facealpha', 0.3);
end

% plot them together, this should form a consistent set
ft_plot_mesh(dfs{1}, 'edgecolor', 'none', 'facecolor', [1 1 1]/2, 'facealpha', 0.3);
ft_plot_mesh(mesh{3}, 'edgecolor', 'none', 'facecolor', [1 1 1]/2, 'facealpha', 0.3);
ft_plot_mesh(surf{3}, 'edgecolor', 'none', 'facecolor', [1 1 1]/2, 'facealpha', 0.3);

dfsall = cellfun(@ft_read_headshape, dfsfile);

