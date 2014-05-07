function test_bug2563

% TEST test_bug2563
% TEST ft_selectdata getdimord

% WALLTIME 00:10:00
% MEM 1gb


nchan = 3;
nfreq = 4;
ntime = 5;

freq = [];
freq.MIspctrm = randn(nchan,nfreq,ntime);
freq.freq = 1:nfreq;
freq.time = 1:ntime;
for i=1:nchan
  freq.label{i} = num2str(i);
end

cfg = [];
cfg.avgoverfreq = 'yes'
output = ft_selectdata(cfg, freq);

assert(isequal(size(output.MIspctrm), [3 1 5]));

disp(freq)
disp(output)



