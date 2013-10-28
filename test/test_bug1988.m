function test_bug1988

% MEM 2gb
% WALLTIME 00:12:37

% TEST test_bug1988 ft_volumesegment ft_prepare_headmodel

%% segmentedmri.mat
% from current version may not match what is on the ftp for tutorials

% as it's called in the BF tutorial
mri=ft_read_mri('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/Subject01.mri');
cfg           = [];
cfg.coordsys  = 'ctf';
segmentedmri_bf  = ft_volumesegment(cfg, mri);

load /home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat
segbf=segmentedmri;clear segmentedmri

% Note to developer: if these assertss fail, is current code wrong, or should
% tutorial be updated?
assert(isequalwithequalnans(segmentedmri_bf.gray,segbf.gray))
assert(isequalwithequalnans(segmentedmri_bf,segbf))

% headmodel_meg, as it's called in the headmodel_meg tutorial
mri=ft_read_mri('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/Subject01.mri');
cfg           = [];
cfg.coordsys  = 'ctf';
cfg.output    = 'brain';
segmentedmri_hm  = ft_volumesegment(cfg, mri);

load /home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_meg/segmentedmri.mat
seghm=segmentedmri;clear segmentedmri

assert(isequalwithequalnans(segmentedmri_hm.brain,seghm.brain))
% length(find(segmentedmri_hm.brain-seghm.brain))
assert(isequalwithequalnans(segmentedmri_hm,seghm))

% headmodel_eeg
load /home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_eeg/segmentedmri.mat
seghe=segmentedmri;clear segmentedmri
mri=ft_read_mri('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/Subject01.mri');
cfg           = [];
cfg.coordsys  = 'ctf';
cfg.output    = {'brain' 'skull' 'scalp'};
segmentedmri_he  = ft_volumesegment(cfg, mri);

assert(isequalwithequalnans(segmentedmri_he.brain,seghe.brain))
% length(find(segmentedmri_he.brain-seghe.brain))
assert(isequalwithequalnans(segmentedmri_he,seghe))



%% vol.mat 
% from current version may not match what is on the ftp for tutorials

load /home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/vol.mat  
% has a vol of type 'nolte' in 'cm' and as an example pnt:
volbf=vol;clear vol
volbf=rmfield(volbf,'cfg');

cfg=[];
cfg.method = 'singleshell';
volbf_new = ft_prepare_headmodel(cfg,segmentedmri_bf);
volbf_new=rmfield(volbf_new,'cfg');

assert(isequalwithequalnans(volbf_new,volbf))


load /home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_meg/vol.mat
% has a vol of type 'singleshell' in 'cm' and as an example pnt:
volhm=vol;clear vol
volhm=rmfield(volhm,'cfg');

cfg=[];
cfg.method = 'singleshell';
volhm_new = ft_prepare_headmodel(cfg,segmentedmri_hm);
volhm_new=rmfield(volhm_new,'cfg');

assert(isequalwithequalnans(volhm_new,volhm))




