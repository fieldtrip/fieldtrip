function test_bug1254

% WALLTIME 00:03:01

% the bug has not been fixed yet, so there is no point in automatically
% executing this test
return

% TEST test_bug1254
% TEST qsubcompile qsubcellfun

% list the possible locations for the data
dataset = {
  '/Users/robert/data/Subject01/Subject01.ds'
  '/home/common/matlab/fieldtrip/data/Subject01.ds'
  };

% pick the available dataset
for i=1:length(dataset)
  if exist(dataset{i}, 'dir')
    dataset = dataset{i};
    break
  end
end

if ~ischar(dataset)
  error('could not locate the test dataset');
end

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


