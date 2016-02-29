function test_bug3075

% MEM 1gb
% WALLTIME 00:20:00

% TEST ft_preprocessing test_bug3075

% the reported issue is that cfg.inputfile and cfg.outputfile do not work

% create some data
data = [];
data.trial = {randn(2,1000)};
data.time  = {(0:999)./1000};
data.label = {'chan01';'chan02'};
data.cfg = [];

filename1 = [tempname,'.mat'];
filename2 = [tempname,'.mat'];

save(filename1, 'data');

% check
cfg = [];
cfg.inputfile = filename1;
%cfg.lpfilter  = 'yes';
%cfg.lpfreq    = 10;
data1 = ft_preprocessing(cfg);

cfg.outputfile = filename2;
ft_preprocessing(cfg);

data2 = ft_preprocessing(cfg);

file1 = load(filename1);
file2 = load(filename2);

assert(isequal(rmfield(data,  'cfg'), rmfield(file1.data, 'cfg'))); % the callinfo may be different
assert(isequal(rmfield(data1, 'cfg'), rmfield(file2.data, 'cfg'))); % the callinfo may be different
assert(isequal(rmfield(data2, 'cfg'), rmfield(file2.data, 'cfg'))); % the callinfo may be different


