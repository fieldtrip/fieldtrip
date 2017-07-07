function inspect_bug1230

% MEM 1500mb
% WALLTIME 00:10:00

% this is not really pertaining to a bug, but to a feature request.

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'));

cfg = [];
cfg.method = 'summary';
cfg.layout = 'CTF151.lay';
data = ft_rejectvisual(cfg, data);
