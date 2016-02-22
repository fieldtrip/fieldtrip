function test_bug3075

% MEM 2gb
% WALLTIME 00:20:00

% TEST ft_preprocessing test_bug3075

% the reported issue is that cfg.inputfile and cfg.outputfile does not work

% create some data
data = [];
data.trial = {randn(2,1000)};
data.time  = {(0:999)./1000};
data.label = {'chan01';'chan02'};

filename1 = [tempname,'.mat'];
filename2 = [tempname,'.mat'];
save(filename1, 'data');

% check
cfg = [];
cfg.inputfile = filename1;
cfg.lpfilter  = 'yes';
cfg.lpfreq    = 10;
datanew = ft_preprocessing(cfg);

cfg.outputfile = filename2;
ft_preprocessing(cfg);


datanew2 = ft_preprocessing(cfg);
tmp = load(filename2);
assert(isequal(datanew2, tmp.data));



