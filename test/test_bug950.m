function test_bug950

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_megrealign test_bug950
% DATA private

% the issue explored here is a reputed crash in megrealign due to a problem
% in the channelposition function.
% in addition balancing is not properly taken into account when creating
% the headmodel for the inverse/forward steps

% load in some data
load(dccnpath('/project/3031000.02/test/latest/raw/meg/preproc_ctf151.mat'));

cfg = [];
cfg.gradient = 'G3BR';
data = ft_denoise_synthetic(cfg, data);

template = data.grad;
template.chanpos(:,3) = template.chanpos(:,3)+1;
template.coilpos(:,3) = template.coilpos(:,3)+1;

cfg = [];
cfg.template{1}    = template;
cfg.inwardshift    = 1;
cfg.headmodel.o    = [0 0 4];
cfg.headmodel.r    = 8;
cfg.headmodel.unit = 'cm';
data2 = ft_megrealign(cfg, data);
