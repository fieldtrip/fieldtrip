function test_bug1297

% MEM 3gb
% WALLTIME 00:30:00
% DEPENDENCY ft_volumesegment

mri_nom = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1297/orig-nomask.mgz'));

mri_nom.coordsys = 'spm'; % this can also be determined with ft_determine_coordsys

cfg           = [];
% cfg.coordsys  = 'spm'; % not supported any more, should be specified in the input data
cfg.output    = {'skullstrip' 'brain'};
seg2           = ft_volumesegment(cfg, mri_nom);
