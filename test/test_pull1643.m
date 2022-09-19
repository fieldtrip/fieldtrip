function test_pull1643

% MEM 8gb
% WALLTIME 00:40:00
% DEPENDENCY ft_connectivity_plm ft_connectivityanalysis

% create some dummy data
ntrl = 100;
nsmp = 1000;
fs   = 1000;
nchan = 3;
for k = 1:ntrl
  trial{1,k} = randn(nchan,nsmp);
  time{1,k}  = (0:999)/fs;
end
for k = 1:nchan
  label{k,1} = sprintf('chan%02d',k);
end
data.trial = trial;
data.time  = time;
data.label = label;

cfg = [];
cfg.method = 'plm';
plm = ft_connectivityanalysis(cfg,data);

% So far the test code only 'checks' whether it runs through, no checks are
% done on the validity of the output structure (in terms of how it looks
% w.r.t. fields etc., nor on the validity of the numeric results)


