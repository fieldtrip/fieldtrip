function test_bug1956

% MEM 1500mb
% WALLTIME 00:60:00

% TEST test_bug1956 
% TEST ft_prepare_sourcemodel volumesmooth

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/nifti/single_subj_T1.nii'));

% up to ~1 November 2013 it was allowed to specify cfg.coordsys, which would then
% get copied over into the input data in ft_sourceplot, ft_volumennormalise and
% ft_volumesegment. From November 2013 onwards, ft_checkdata is used with
% hascoordsys=yes.

mri.coordsys = 'ctf'; % this can also be determined with ft_determine_coordsys

cfg=[];
% cfg.coordsys = 'ctf'; % not supported any more, should be specified in the input data
tpm = ft_volumesegment(cfg,mri);  % gray, white, csf tissue prob. map

cfg=[];
cfg.mri=tpm;
grid01 = ft_prepare_sourcemodel(cfg,tpm);

cfg=[];
cfg.mri=tpm;
grid02 = ft_prepare_sourcemodel(cfg,tpm);

cfg=[];
cfg.mri=tpm;
grid03 = ft_prepare_sourcemodel(cfg,tpm);

grid01 = rmfield(grid01, 'cfg');
grid02 = rmfield(grid02, 'cfg');
grid03 = rmfield(grid03, 'cfg');

% make sure that they are the same, and not super-smoohted
assert(isequal(grid01,grid02))
assert(isequal(grid01,grid03))
