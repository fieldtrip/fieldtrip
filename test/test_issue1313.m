function test_issue1313

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY resampledesign

%%

[ftver, ftpath] = ft_version;
cd(fullfile(ftpath, 'private'));

%%

design = [
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  ];

cfg = [];
cfg.uvar = 1;
cfg.ivar = 2;
cfg.resampling = 'permutation';
cfg.numrandomization = 'all';
[resample] = resampledesign(cfg, design);

assert(size(resample,1)==1024);

%%

ft_warning('nothing')

design = [
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  ];

cfg = [];
cfg.uvar = 1;
cfg.ivar = 2;
cfg.resampling = 'permutation';
cfg.numrandomization = 600;
[resample] = resampledesign(cfg, design);

assert(size(resample,1)==600);

% ensure that the requested warning is given
w = ft_warning('last');
assert(startsWith(w.message, 'the number of randomizations'));

%%
% 10 subjects, 2 conditions

ft_warning('nothing')

design = [
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  ];

cfg = [];
cfg.uvar = 1;
cfg.ivar = 2;
cfg.resampling = 'permutation';
cfg.numrandomization = 1500;
[resample] = resampledesign(cfg, design);

assert(size(resample,1)==1500);

% ensure that the requested warning is given
w = ft_warning('last');
assert(startsWith(w.message, 'the number of randomizations'));
