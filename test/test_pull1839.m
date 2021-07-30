function test_pull1839

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_timelockbaseline

%%

% create some data
data   = [];
nsmp   = [100 100 100 100];
offset = [-10 -10 -10 -10];
for k = 1:numel(nsmp)
  data.trial{1,k} = randn(2,nsmp(k));
  data.time{1,k}  = ((0:(nsmp(k)-1))+offset(k))./100;
end
data.label = {'chan01';'chan02'};

% it should also work for component data
cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);

% it should also work for timelock data
cfg = [];
cfg.keeptrials = 'yes';
tlck = ft_timelockanalysis(cfg, data);

cfg1 = [];
cfg1.demean = 'yes';
cfg1.baselinewindow = [-0.1 0];
dataout1 = ft_preprocessing(cfg1, data);
compout1 = ft_preprocessing(cfg1, comp);
tlckout1 = ft_preprocessing(cfg1, tlck);

cfg2 = [];
cfg2.baseline = [-0.1 0];
dataout2 = ft_timelockbaseline(cfg2, data);
compout2 = ft_timelockbaseline(cfg2, comp);
tlckout2 = ft_timelockbaseline(cfg2, tlck);

% this should work just fine
[ok,message] = isalmostequal(rmfield(dataout1,'cfg'),rmfield(dataout2,'cfg'),'abstol',100*eps);
assert(ok);
[ok,message] = isalmostequal(rmfield(compout1,'cfg'),removefields(compout2,{'cfg' 'topodimord' 'unmixingdimord'}),'abstol',100*eps);
assert(ok);
[ok,message] = isalmostequal(rmfield(tlckout1,'cfg'),rmfield(tlckout2,'cfg'),'abstol',100*eps);
assert(ok);

% create some other data
data   = [];
nsmp   = [100 110 120 130];
offset = [-10 -10 -10 -10];
for k = 1:numel(nsmp)
  data.trial{1,k} = randn(2,nsmp(k));
  data.time{1,k}  = ((0:(nsmp(k)-1))+offset(k))./100;
end
data.label = {'chan01';'chan02'};

% it should also work for component data
cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);

% it should also work for timelock data
cfg = [];
cfg.keeptrials = 'yes';
tlck = ft_timelockanalysis(cfg, data);

cfg1 = [];
cfg1.demean = 'yes';
cfg1.baselinewindow = [-0.1 0];
dataout1 = ft_preprocessing(cfg1, data);
compout1 = ft_preprocessing(cfg1, comp);
tlckout1 = ft_preprocessing(cfg1, tlck);

cfg2 = [];
cfg2.baseline = [-0.1 0];
dataout2 = ft_timelockbaseline(cfg2, data);
compout2 = ft_timelockbaseline(cfg2, comp);
tlckout2 = ft_timelockbaseline(cfg2, tlck);

% this inserts nans for the shorter trials in the converted data objects,
% dataout2/compout2
[ok,message] = isalmostequal(dataout1.trial{1},dataout2.trial{1}(:,1:100),'abstol',100*eps);
assert(ok);
assert(numel(dataout1.time{1}) ~= numel(dataout2.time{1}));
[ok,message] = isalmostequal(compout1.trial{1},compout2.trial{1}(:,1:100),'abstol',100*eps);
assert(ok);
[ok,message] = isalmostequal(removefields(tlckout1,{'cfg' 'sampleinfo'}),rmfield(tlckout2,'cfg'),'abstol',100*eps);
assert(ok);

% create yet some other data, now including a trial that does not have data
% for the requested basline window
data   = [];
nsmp   = [100 100 100 100];
offset = [-10 -10 -10 11];
for k = 1:numel(nsmp)
  data.trial{1,k} = randn(2,nsmp(k)); % add offset in amplitude values to be able to evaluate later
  if k<4,data.trial{1,k}(:,1:10) = 4; end
  data.time{1,k}  = ((0:(nsmp(k)-1))+offset(k))./100;
  
end
data.label = {'chan01';'chan02'};

% it should also work for component data
cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);

% it should also work for timelock data
cfg = [];
cfg.keeptrials = 'yes';
tlck = ft_timelockanalysis(cfg, data);

cfg1 = [];
cfg1.demean = 'yes';
cfg1.baselinewindow = [-0.1 0];
dataout1 = ft_preprocessing(cfg1, data);
compout1 = ft_preprocessing(cfg1, comp);
tlckout1 = ft_preprocessing(cfg1, tlck);

cfg2 = [];
cfg2.baseline = [-0.1 0];
dataout2 = ft_timelockbaseline(cfg2, data);
compout2 = ft_timelockbaseline(cfg2, comp);
tlckout2 = ft_timelockbaseline(cfg2, tlck);

% this inserts nans for the shorter trials in the converted data objects,
% dataout2/compout2
[ok,message] = isalmostequal(dataout1.trial{4},dataout2.trial{4}(:,22:121),'abstol',100*eps);
assert(ok);
assert(numel(dataout1.time{1}) ~= numel(dataout2.time{1}));
[ok,message] = isalmostequal(compout1.trial{4},compout2.trial{4}(:,22:121),'abstol',100*eps);
assert(ok);

