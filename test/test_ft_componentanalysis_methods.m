function test_ft_componentanalysis_methods

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_componentanalysis, ft_rejectcomponent, ft_checkdata
% TEST runica, fastica, sobi, jader, dss, parafac

% note 3) binica does not run on all computers, which is why it is disabled
% note 9) parafac depends on an external toolbox, which I don't have access to at the moment 

% construct a simple and minimalistic raw data representation
nchan  = 10;
ntime  = 1000;
ntrial = 5;
data = [];
data.fsample = 1;
for i=1:nchan
  data.label{i} = sprintf('chan%d', i);
end
for i=1:ntrial
  data.trial{i} = sqrt(abs(randn(nchan,ntime)));
  data.time{i} = 1:ntime;
end

% run the decomposition with various methods
cfg = [];
cfg.method = 'fastica';
comp1f = ft_componentanalysis(cfg, data);
cfg.numcomponent = nchan-1;
comp1s = ft_componentanalysis(cfg, data);

cfg = [];
cfg.method = 'runica';
cfg.runica.maxsteps = 50;
comp2f = ft_componentanalysis(cfg, data);
cfg.numcomponent = nchan-1;
comp2s = ft_componentanalysis(cfg, data);

% cfg = [];
% cfg.method = 'binica';
% comp3f = ft_componentanalysis(cfg, data);
% cfg.numcomponent = nchan-1;
% comp3s = ft_componentanalysis(cfg, data);

cfg = [];
cfg.method = 'jader';
comp4f = ft_componentanalysis(cfg, data);
cfg.numcomponent = nchan-1;
comp4s = ft_componentanalysis(cfg, data);

cfg = [];
cfg.method = 'varimax';
comp5f = ft_componentanalysis(cfg, data);
cfg.numcomponent = nchan-1;
comp5s = ft_componentanalysis(cfg, data);

cfg = [];
cfg.method = 'cca';
comp6f = ft_componentanalysis(cfg, data);
cfg.numcomponent = nchan-1;
comp6s = ft_componentanalysis(cfg, data);

cfg = [];
cfg.method = 'pca';
comp7f = ft_componentanalysis(cfg, data);
cfg.numcomponent = nchan-1;
comp7s = ft_componentanalysis(cfg, data);

cfg = [];
cfg.method = 'svd';
comp8f = ft_componentanalysis(cfg, data);
cfg.numcomponent = nchan-1;
comp8s = ft_componentanalysis(cfg, data);

% cfg = [];
% cfg.method = 'parafac';
% comp9f = ft_componentanalysis(cfg, data);
% cfg.numcomponent = nchan-1;
% comp9s = ft_componentanalysis(cfg, data);

cfg = [];
cfg.method = 'dss';
comp10f = ft_componentanalysis(cfg, data);
cfg.numcomponent = nchan-1;
comp10s = ft_componentanalysis(cfg, data);

cfg = [];
cfg.method = 'sobi';
comp11f = ft_componentanalysis(cfg, data);
cfg.numcomponent = nchan-1;
comp11s = ft_componentanalysis(cfg, data);

% reproject the data to the channel level
cfg = [];
[dataout] = ft_rejectcomponent(cfg, comp1f);
[dataout] = ft_rejectcomponent(cfg, comp2f);
% [dataout] = ft_rejectcomponent(cfg, comp3f);
[dataout] = ft_rejectcomponent(cfg, comp4f);
[dataout] = ft_rejectcomponent(cfg, comp5f);
[dataout] = ft_rejectcomponent(cfg, comp6f);
[dataout] = ft_rejectcomponent(cfg, comp7f);
[dataout] = ft_rejectcomponent(cfg, comp8f);
% [dataout] = ft_rejectcomponent(cfg, comp9f);
[dataout] = ft_rejectcomponent(cfg, comp10f);
[dataout] = ft_rejectcomponent(cfg, comp11f);

cfg = [];
[dataout] = ft_rejectcomponent(cfg, comp1s);
[dataout] = ft_rejectcomponent(cfg, comp2s);
% [dataout] = ft_rejectcomponent(cfg, comp3s);
[dataout] = ft_rejectcomponent(cfg, comp4s);
[dataout] = ft_rejectcomponent(cfg, comp5s);
[dataout] = ft_rejectcomponent(cfg, comp6s);
[dataout] = ft_rejectcomponent(cfg, comp7s);
[dataout] = ft_rejectcomponent(cfg, comp8s);
% [dataout] = ft_rejectcomponent(cfg, comp9s);
[dataout] = ft_rejectcomponent(cfg, comp10s);
[dataout] = ft_rejectcomponent(cfg, comp11s);

cfg = [];
cfg.component = [1 3];
[dataout] = ft_rejectcomponent(cfg, comp1f);
[dataout] = ft_rejectcomponent(cfg, comp2f);
% [dataout] = ft_rejectcomponent(cfg, comp3f);
[dataout] = ft_rejectcomponent(cfg, comp4f);
[dataout] = ft_rejectcomponent(cfg, comp5f);
[dataout] = ft_rejectcomponent(cfg, comp6f);
[dataout] = ft_rejectcomponent(cfg, comp7f);
[dataout] = ft_rejectcomponent(cfg, comp8f);
% [dataout] = ft_rejectcomponent(cfg, comp9f);
[dataout] = ft_rejectcomponent(cfg, comp10f);
[dataout] = ft_rejectcomponent(cfg, comp11f);

cfg = [];
cfg.component = [1 3];
[dataout] = ft_rejectcomponent(cfg, comp1s);
[dataout] = ft_rejectcomponent(cfg, comp2s);
% [dataout] = ft_rejectcomponent(cfg, comp3s);
[dataout] = ft_rejectcomponent(cfg, comp4s);
[dataout] = ft_rejectcomponent(cfg, comp5s);
[dataout] = ft_rejectcomponent(cfg, comp6s);
[dataout] = ft_rejectcomponent(cfg, comp7s);
[dataout] = ft_rejectcomponent(cfg, comp8s);
% [dataout] = ft_rejectcomponent(cfg, comp9s);
[dataout] = ft_rejectcomponent(cfg, comp10s);
[dataout] = ft_rejectcomponent(cfg, comp11s);

% perform some checks on the decomposed data properties
assert(length(comp1f.label)==nchan);
assert(length(comp2f.label)==nchan);
% assert(length(comp3f.label)==nchan);
assert(length(comp4f.label)==nchan);
assert(length(comp5f.label)==nchan);
assert(length(comp6f.label)==nchan);
assert(length(comp7f.label)==nchan);
assert(length(comp8f.label)==nchan);
assert(length(comp8f.label)==nchan);
% assert(length(comp9f.label)==nchan);
assert(length(comp10f.label)==nchan);
assert(length(comp11f.label)==nchan);

% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=903
assert(length(comp1s.label)==nchan-1);
assert(length(comp2s.label)==nchan-1);
% assert(length(comp3s.label)==nchan);
assert(length(comp4s.label)==nchan-1);
assert(length(comp5s.label)==nchan-1);
assert(length(comp6s.label)==nchan-1);
assert(length(comp7s.label)==nchan-1);
assert(length(comp8s.label)==nchan-1);
% assert(length(comp9s.label)==nchan-1);
assert(length(comp10s.label)==nchan-1);
assert(length(comp11s.label)==nchan-1);



