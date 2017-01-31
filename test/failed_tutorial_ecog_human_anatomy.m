function failed_tutorial_ecog_human_anatomy

% MEM 5gb
% WALLTIME 00:10:00

% TEST test_tutorial_ecog_human_anatomy
% TEST

datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/human_ecog');
subj = 'SubjectUCI29';
subj_dir = fullfile(datadir, subj);
cd(subj_dir)

% use matlab functions that come with freesurfer
% the subset in FieldTrip is not enough. The problem is freesurfer's
% license (does it allow to copy their functions to GPLv3 toolbox?)
[~, recon_file] = system('which recon-all');
freesurfer_dir = fileparts(fileparts(recon_file));

addpath(fullfile(freesurfer_dir, 'matlab'))

load('SubjectUCI29_elec.mat')

% TUTORIAL: Surface rendering
mesh = ft_read_headshape('SubjectUCI29_lh.pial');
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none')
view([-90 25])
lighting gouraud
material shiny
camlight

% plot electrodes
hs = ft_plot_sens(elec, 'style', 'ko');
set(hs, 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% smooth outer surface
cfg = [];
cfg.headshape = 'SubjectUCI29_lh.pial';
cfg.method = 'cortexhull';
smooth_mesh = ft_prepare_mesh(cfg);  % this function uses Freesurfer!

% plot smoothed mesh
ft_plot_mesh(smooth_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none')
view([-90 25])
lighting gouraud
material shiny
camlight

% snap to pial

cfg = [];
cfg.method = 'headshape';
cfg.warp = 'dykstra2012';
cfg.channel = {'POL *G*'};
cfg.elec = elec;
cfg.headshape = smooth_mesh;
cfg.feedback = 'yes';

% snap to pial
elec_snapped = ft_electroderealign(cfg);

% plot snapped elec
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none')
view([-90 25])
lighting gouraud
material shiny
camlight

% plot electrodes
hs = ft_plot_sens(elec_snapped, 'style', 'ko');
set(hs, 'MarkerFaceColor', 'k', 'MarkerSize', 6);
