function test_bug1254

% MEM 1500mb
% WALLTIME 00:10:00

% the bug has not been fixed yet, so there is no point in automatically
% executing this test
return

% TEST qsubcompile qsubcellfun

dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');

cfg = [];
cfg.dataset = dataset;
cfg.trl = [
     1  900 0
   901 1800 0
  1801 2700 0
  ];
data1 = ft_preprocessing(cfg);

% this is where it gets critical
compiledfun = qsubcompile('ft_preprocessing');
data2 = qsubcellfun(compiledfun, {cfg}, 'memreq', 1024^3, 'timreq', 300);

% this also failed
compiledfun = qsubcompile('istrue');
value = qsubcellfun(compiledfun, {'yes'}, 'memreq', 1024^3, 'timreq', 300);


