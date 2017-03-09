function test_bug2563

% TEST ft_selectdata getdimord

% WALLTIME 00:10:00
% MEM 1gb


nsubj = 3;
nchan = 4;
nfreq = 5;

freq = [];
freq.dimord = 'subj_chan_freq';
freq.MIspctrm = randn(nsubj,nchan,nfreq);
freq.freq = 1:nfreq;
for i=1:nchan
  freq.label{i} = num2str(i);
end

cfg = [];
cfg.avgoverfreq = 'yes';
output = ft_selectdata(cfg, freq);

assert(isequal(size(output.MIspctrm), [3 4]));

disp(freq)
disp(output)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following was reported as http://bugzilla.fcdonders.nl/show_bug.cgi?id=2561
% but seems to be related, i.e. both have subjects in the data

timelock = [];
timelock.individual = randn(7,204,1500);
timelock.time = 1:1500;
for i=1:204
  timelock.label{i} = num2str(i);
end
timelock.dimord = 'subj_chan_time';

cfg = [];
output = ft_selectdata(cfg, timelock);

