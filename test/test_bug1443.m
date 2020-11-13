function test_bug1443

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_rejectcomponent ft_componentanalysis

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'));

cfg = [];
cfg.method = 'fastica';
cfg.numcomponent = 14; % to make it go fast
cfg.randomseed = 13; % so we get the same output each time
comp = ft_componentanalysis(cfg,data);

assert(ft_senstype(comp,'ctf151'))

cfg = [];
cfg.component = 2; % chosen randomly
rej1 = ft_rejectcomponent(cfg, comp, data);
rej2 = ft_rejectcomponent(cfg, comp);

norm(rej2.grad.tra-rej1.grad.tra)/norm(rej2.grad.tra);
figure; imagesc(rej2.grad.tra - rej1.grad.tra); caxis([-1 1])

load standard_sourcemodel3d10mm
load standard_singleshell

cfg = [];
cfg.sourcemodel = sourcemodel;
cfg.headmodel = vol;

cfg.grad = rej1.grad;
grid1 = ft_prepare_leadfield(cfg, rej1);

cfg.grad = rej2.grad;
grid2 = ft_prepare_leadfield(cfg, rej2);

assert(~isequaln(grid1.leadfield,grid2.leadfield))
