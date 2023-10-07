function test_issue2021

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_preprocessing
% DATA no

data  = randn(10);
fname = [tempname '.mat'];
save(fname, 'data');

cfg         = [];
cfg.dataset = fname;
try
  data        = ft_preprocessing(cfg);
catch
  assert(contains(lasterr, 'This function requires raw data as input, see ft_datatype_raw.'));
end
