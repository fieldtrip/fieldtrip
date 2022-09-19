function test_bug3403

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

ntime   = 1000;
ntrial  = 10;
nchan   = 3;
fsample = 1000;

data = [];
for i=1:nchan
  data.label{i} = num2str(i);
end

for i=1:ntrial
  data.time{i} = (1:ntime)/fsample;
  data.trial{i} = randn(nchan, ntime);
end

% insert an artifact
data.trial{2}(2,:) = linspace(-20,20,ntime);

%%

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.range = 10;
[cfg1, artefact1] = ft_artifact_threshold(cfg, data);

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.min = -10;
[cfg2, artefact2] = ft_artifact_threshold(cfg, data);

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.max = 10;
[cfg3, artefact3] = ft_artifact_threshold(cfg, data);

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.min = -10;
cfg.artfctdef.threshold.max =  10;
[cfg4, artefact4] = ft_artifact_threshold(cfg, data);

%%

cfg1.artfctdef.reject = 'none';
data1none = ft_rejectartifact(cfg1, data);
assert(isequal(data1none.trial, data.trial));

cfg1.artfctdef.reject = 'partial';
data1partial = ft_rejectartifact(cfg1, data);
assert(numel(data1partial.trial)~=ntrial);

cfg1.artfctdef.reject = 'nan';
data1nan = ft_rejectartifact(cfg1, data);
assert(numel(data1nan.trial)==ntrial);

cfg1.artfctdef.reject = 'complete';
data1complete = ft_rejectartifact(cfg1, data);
assert(numel(data1complete.trial)~=ntrial);

%%

cfg2.artfctdef.reject = 'none';
data2none = ft_rejectartifact(cfg2, data);
assert(isequal(data2none.trial, data.trial));

cfg2.artfctdef.reject = 'partial';
data2partial = ft_rejectartifact(cfg2, data);
assert(numel(data2partial.trial)==ntrial);

cfg2.artfctdef.reject = 'nan';
data2nan = ft_rejectartifact(cfg2, data);
assert(numel(data2nan.trial)==ntrial);

cfg2.artfctdef.reject = 'complete';
data2complete = ft_rejectartifact(cfg2, data);
assert(numel(data2complete.trial)~=ntrial);

%%

cfg3.artfctdef.reject = 'none';
data3none = ft_rejectartifact(cfg3, data);
assert(isequal(data3none.trial, data.trial));

cfg3.artfctdef.reject = 'partial';
data3partial = ft_rejectartifact(cfg3, data);
assert(numel(data3partial.trial)==ntrial);

cfg3.artfctdef.reject = 'nan';
data3nan = ft_rejectartifact(cfg3, data);
assert(numel(data3nan.trial)==ntrial);

cfg3.artfctdef.reject = 'complete';
data3complete = ft_rejectartifact(cfg3, data);
assert(numel(data3complete.trial)~=ntrial);

%%

cfg4.artfctdef.reject = 'none';
data4none = ft_rejectartifact(cfg4, data);
assert(isequal(data4none.trial, data.trial));

cfg4.artfctdef.reject = 'partial';
data4partial = ft_rejectartifact(cfg4, data);
assert(numel(data4partial.trial)==ntrial);

cfg4.artfctdef.reject = 'nan';
data4nan = ft_rejectartifact(cfg4, data);
assert(numel(data4nan.trial)==ntrial);

cfg4.artfctdef.reject = 'complete';
data4complete = ft_rejectartifact(cfg4, data);
assert(numel(data4complete.trial)~=ntrial);
