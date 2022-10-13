function test_ft_freqanalysis

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqanalysis

ntrl = 1000;
nchan = 32;
fs = 500;
start_time = -1; % seconds
end_time = 2.5; % seconds
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.label = cellstr(num2str((1:nchan).'));
for i=1:ntrl
  data.time{i} = linspace(start_time, end_time, nsamples);
  data.trial{i} = randn(nchan,nsamples);
end

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [4 30];
freq = ft_freqanalysis(cfg, data);
