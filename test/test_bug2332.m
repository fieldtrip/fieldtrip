function test_bug2332

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug2332
% TEST ft_apply_montage ft_componentanalysis ft_rejectcomponent

%% read some data

pwdir = pwd;

cd('/home/common/matlab/fieldtrip/data/');
cfg = [];
cfg.dataset = 'Subject01.ds';
cfg.trl     = [1 1200 0];
cfg.continuous = 'yes';
cfg.channel = 'MEG';
data = ft_preprocessing(cfg);

%% do componentanalysis

cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);

%% verify grad has been updated successfully

assert(isequal(comp.grad.chanposorg, data.grad.chanpos));
assert(all(isnan(comp.grad.chanpos(:))));

%% rejectcomponent

cfg = [];
cfg.component = 2;
reject = ft_rejectcomponent(cfg, comp);

%% megplanar

cfg = [];
cfg.method = 'template';
cfg.template = 'ctf151_neighb.mat';
neighbours = ft_prepare_neighbours(cfg);

cfg = [];
cfg.neighbours = neighbours;
cfg.method = 'sincos';
datplan = ft_megplanar(cfg, reject);

%% verify again that grad was updated successfully

assert(all(isnan(datplan.grad.chanpos(:))));
assert(~any(isnan(datplan.grad.chanposorg(:))));

% *org fields should still reflect the original physical positions
assert(isequal(datplan.grad.chanposorg, data.grad.chanpos));

%% simple check that combineplanar still works

datcmb = ft_combineplanar([], datplan);

% now chanpos should be nans
assert(all(isnan(datcmb.grad.chanpos(:))));
% FIXME do we want a chanposorg on combined planar data as well?

%% move back to original working dir

cd(pwdir);

