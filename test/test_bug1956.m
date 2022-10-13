function test_bug1956

% WALLTIME 02:00:00
% MEM 3gb
% DEPENDENCY ft_volumesegment ft_prepare_sourcemodel volumesmooth

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.mri'));
mri.coordsys = 'ctf'; % this can also be determined with ft_determine_coordsys

cfg = [];
tpm = ft_volumesegment(cfg,mri);  % gray, white, csf tissue prob. map

% there are some spm functions which change their own input

cfg = [];
cfg.mri = tpm;
grid01 = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.mri = tpm;
grid02 = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.mri = tpm;
grid03 = ft_prepare_sourcemodel(cfg);

grid01 = rmfield(grid01, 'cfg');
grid02 = rmfield(grid02, 'cfg');
grid03 = rmfield(grid03, 'cfg');

% make sure that they are the same, and not super-smoohted
assert(isequal(grid01,grid02))
assert(isequal(grid01,grid03))
