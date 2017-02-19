function failed_bug1988

% MEM 3gb
% WALLTIME 00:45:00

% TEST test_bug1988 ft_volumesegment ft_prepare_headmodel

%% segmentedmri.mat
% from current version may not match what is on the ftp for tutorials

% as it's called in the BF tutorial
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/Subject01.mri'));

mri.coordsys = 'ctf'; % this can also be determined with ft_determine_coordsys

cfg = [];
% cfg.coordsys = 'ctf'; % not supported any more, should be specified in the input data
segmentedmri_bf = ft_volumesegment(cfg, mri);

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat'));
reference_bf = segmentedmri;clear segmentedmri

% Note to developer: if these assertss fail, is current code wrong, or should
% tutorial be updated?
assert(isequaln(segmentedmri_bf.gray,reference_bf.gray))
assert(isequaln(rmfield(segmentedmri_bf, 'cfg'),rmfield(reference_bf, 'cfg')))

% headmodel_meg, as it's called in the headmodel_meg tutorial
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/Subject01.mri'));

mri.coordsys = 'ctf'; % this can also be determined with ft_determine_coordsys

cfg = [];
% cfg.coordsys = 'ctf'; % not supported any more, should be specified in the input data
cfg.output = 'brain';
segmentedmri_hm = ft_volumesegment(cfg, mri);

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_meg/segmentedmri.mat'));
reference_hm = segmentedmri;clear segmentedmri

reference_hm = tryrmfield(reference_hm, 'cfg');
segmentedmri_hm = tryrmfield(segmentedmri_hm, 'cfg');

assert(isequaln(segmentedmri_hm.brain,reference_hm.brain))
assert(isequaln(segmentedmri_hm,reference_hm))

% headmodel_eeg
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_eeg/segmentedmri.mat'));
reference_he = segmentedmri;clear segmentedmri

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/Subject01.mri'));

mri.coordsys = 'ctf'; % this can also be determined with ft_determine_coordsys

cfg = [];
% cfg.coordsys = 'ctf'; % not supported any more, should be specified in the input data
cfg.output = {'brain' 'skull' 'scalp'};
segmentedmri_he = ft_volumesegment(cfg, mri);

reference_he = tryrmfield(reference_he, 'cfg');
segmentedmri_he = tryrmfield(segmentedmri_he, 'cfg');

assert(isequaln(segmentedmri_he.brain,reference_he.brain))
assert(isequaln(segmentedmri_he,reference_he))

%% vol.mat
% from current version may not match what is on the ftp for tutorials

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/vol.mat'));
% has a vol of type 'nolte' in 'cm' and as an example pnt:
volbf = vol; clear vol

cfg = [];
cfg.method = 'singleshell';
volbf_new = ft_prepare_headmodel(cfg,segmentedmri_bf);

volbf = tryrmfield(volbf, 'cfg');
volbf_new = tryrmfield(volbf_new,'cfg');
assert(isequaln(volbf_new,volbf))

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_meg/vol.mat'));
% has a vol of type 'singleshell' in 'cm' and as an example pnt:
volhm = vol; clear vol

cfg = [];
cfg.method = 'singleshell';
volhm_new = ft_prepare_headmodel(cfg,segmentedmri_hm);
volhm_new = rmfield(volhm_new,'cfg');

volhm = rmfield(volhm,'cfg');
volhm_new = rmfield(volhm_new,'cfg');
assert(isequaln(volhm_new,volhm))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = tryrmfield(s, f)
if isfield(s, f)
  s = rmfield(s, f);
end


