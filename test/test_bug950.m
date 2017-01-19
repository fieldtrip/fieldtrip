function test_bug950

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_megrealign test_bug950

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% the issue explored here is a reputed crash in megrealign due to a problem
% in the channelposition function.
% in addition balancing is not properly taken into account when creating
% the headmodel for the inverse/forward steps

% load in some data
load(dccnpath(fullfile('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/','preproc_ctf151')));

cfg = [];
cfg.gradient = 'G3BR';
data = ft_denoise_synthetic(cfg, data);

template = data.grad;
template.chanpos(:,3) = template.chanpos(:,3)+1;
template.coilpos(:,3) = template.coilpos(:,3)+1;

cfg = [];
cfg.template{1} = template;
cfg.inwardshift = 1;
cfg.vol.o    = [0 0 4];
cfg.vol.r    = 8;
cfg.vol.unit = 'cm';
data2 = ft_megrealign(cfg, data);
