function test_bug1230

% MEM 1500mb
% WALLTIME 00:10:00

% this is not really pertaining to a bug, but to a feature request.
% therefore the code is commented out, the more so, because the function that is called
% will enter the interactive mode, and will never end.

%load('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat');
%
%cfg = [];
%cfg.method = 'summary';
%cfg.layout = 'CTF151.lay';
%data = ft_rejectvisual(cfg, data);
