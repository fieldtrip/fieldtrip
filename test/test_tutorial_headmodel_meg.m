function test_tutorial_headmodel_meg(datadir)

% MEM  2000mb
% WALLTIME 00:27:44

% TEST test_tutorial_headmodel_meg
% TEST ft_read_mri ft_volumesegment ft_prepare_headmodel ft_plot_vol
% TEST ft_convert_units ft_read_sens ft_plot_sens 
% intial version by Lilla Magyari

if nargin==0
  datadir = '/home/common/matlab/fieldtrip/data/';
end

mri = ft_read_mri([datadir,'ftp/tutorial/beamformer/Subject01.mri']);

cfg           = [];
cfg.coordsys  = 'ctf';
cfg.output    = 'brain';
segmentedmri  = ft_volumesegment(cfg, mri);

% check if segmentation is equivalent with segmentation on the ftp site

segmentedmri2 = load([datadir,'ftp/tutorial/headmodel_meg/segmentedmri']);
segmentedmri=rmfield(segmentedmri,'cfg');
segmentedmri2=rmfield(segmentedmri2.segmentedmri,'cfg');
assert(isequal(segmentedmri2,segmentedmri),'The segmentation does not match the segmentation stored on the ftp site');

%

cfg = [];
cfg.method='singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);

% check if vol is equivalent with vol on the ftp site

vol2 = load([datadir,'ftp/tutorial/headmodel_meg/vol']);
vol=rmfield(vol,'cfg');
vol2=rmfield(vol2.vol,'cfg');
vol = ft_convert_units(vol,'mm');
assert(identical(vol,vol2,'abstol',0.0001),'The headmodel does not match the headmodel stored on the ftp site.');

% 
sens = ft_read_sens([datadir,'/Subject01.ds']);

vol = ft_convert_units(vol,'cm');

figure
ft_plot_sens(sens, 'style', '*b');

hold on
ft_plot_vol(vol);






