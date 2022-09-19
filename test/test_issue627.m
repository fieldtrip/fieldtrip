function test_issue627

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_preprocessing preproc ft_apply_montage

%%
nchan = 10;
ntrial = 20;
ntime = 1000;
fsample = 1000;

data = [];
for i=1:nchan
  data.label{i} = num2str(i);
end
for i=1:ntrial
  data.trial{i} = randn(nchan, ntime);
  data.time{i} = (1:ntime)/fsample;
end

%%

cfg = [];
data0 = ft_preprocessing(cfg, data);

assert(isequal(data0.label, data.label));
assert(isequal(data0.trial, data.trial));

%%

cfg = [];
cfg.refmethod = 'avg';
cfg.reref = 'no';  % NOTHING SHOULD HAPPEN
cfg.refchannel = 'all';
data1 = ft_preprocessing(cfg, data);

assert(isequal(data1.label, data.label));
assert(isequal(data1.trial, data.trial));

%%

cfg = [];
cfg.refmethod = 'avg';
cfg.reref = 'yes';
cfg.refchannel = 'all';
data2 = ft_preprocessing(cfg, data);

assert( isequal(data2.label, data.label));
assert(~isequal(data2.trial, data.trial));

%%

cfg = [];
cfg.refmethod = 'median';
cfg.reref = 'yes';
cfg.refchannel = 'all';
data3 = ft_preprocessing(cfg, data);

assert( isequal(data3.label, data.label));
assert(~isequal(data3.trial, data.trial));
assert(~isequal(data3.trial, data2.trial)); % they should not be the same

%%

cfg = [];
cfg.refmethod = 'bipolar';
cfg.reref = 'yes';
cfg.refchannel = 'all';
data3 = ft_preprocessing(cfg, data);

assert(~isequal(data3.label, data.label));
assert(~isequal(data3.trial, data.trial));
assert(isequal(data3.trial{1}, data0.trial{1}(1:(nchan-1),:) - data0.trial{1}(2:nchan,:)));
