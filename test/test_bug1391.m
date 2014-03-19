function test_bug1391

% MEM 1500mb
% WALLTIME 00:10:00

% this functions tests some potential issues with numerical accuracy
% of the time axis in preproc, which converts the time axis into a 
% fsample and offset, and then back into a time axis again.

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
cfg.trl     = [10001 11200 -600];
cfg.channel = 'MLC';
cfg.continuous = 'yes';
data1 = ft_preprocessing(cfg);

cfg.dftfilter = 'yes';
cfg.padding   = 10;
cfg.channel   = 'MRC';
cfg.continuous = 'yes';
data2 = ft_preprocessing(cfg);

if ~all(data1.time{1}==data2.time{1})
  error('numerical round off issues detected in the time axes');
end

