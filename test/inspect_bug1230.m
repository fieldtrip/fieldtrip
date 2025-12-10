function inspect_bug1230

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY
% DATA private

% this is not really pertaining to a bug, but to a feature request.

load(dccnpath('/project/3031000.02/test/latest/raw/meg/preproc_ctf151.mat'));

cfg = [];
cfg.method = 'summary';
cfg.layout = 'CTF151.lay';
data = ft_rejectvisual(cfg, data);
