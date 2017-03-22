function test_bug2332

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_apply_montage ft_componentanalysis ft_rejectcomponent

%% read some data

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
cfg.trl     = [1 1200 0];
cfg.continuous = 'yes';
cfg.channel = {'meg', 'megref'};
data = ft_preprocessing(cfg);

%% convert to 3rd order gradients

cfg = [];
cfg.gradient = 'G3BR';
synthetic = ft_denoise_synthetic(cfg, data);

%% do componentanalysis

cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, synthetic);

%% verify grad has been updated successfully

assert(~all(isnan(comp.grad.chanpos(:)))); % as of 8 March 2017 the chanpos is updated

[sel1, sel2] = match_str(comp.grad.labelold, data.grad.label);
assert(isequal(comp.grad.chanposold(sel1,:), data.grad.chanpos(sel2,:))); % old fields should still reflect the original physical positions

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

assert(~all(isnan(datplan.grad.chanpos(:)))); % as of 8 March 2017 the chanpos is updated

[sel1, sel2] = match_str(datplan.grad.labelold, data.grad.label);
assert(isequal(datplan.grad.chanposold(sel1,:), data.grad.chanpos(sel2,:))); % old fields should still reflect the original physical positions

%% simple check that combineplanar still works

cfg = [];
datcmb = ft_combineplanar(cfg, datplan);

% it is not fully clear what the expected channelpositions are at this moment

