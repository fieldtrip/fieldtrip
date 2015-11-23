function test_bug2976

% WALLTIME 00:10:00
% MEM 1gb

%%

fsample = 1000;
ntrial = 20;
nsample = 1000;
nchan = 10;

data = [];
for i=1:nchan
  data.label{i} = sprintf('%d', i);
end
for i=1:ntrial
  data.time{i}  = (1:nsample)./fsample;
  data.trial{i} = randn(nchan, nsample);
  data.trialinfo(i,1) = i;
  data.sampleinfo(i,1) = (i-1)*nsample + 1;
  data.sampleinfo(i,2) = (i  )*nsample;
end

%%

cfg = [];
cfg.artfctdef.manual.artifact = [500 1500];
clean = ft_rejectartifact(cfg, data);

%%

assert(numel(clean.trial)==18);
assert(numel(clean.trialinfo)==18);
assert(isequal(clean.trialinfo(:), (3:20)'));


