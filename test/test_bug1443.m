function test_bug1443

% MEM 1gb
% WALLTIME 00:04:44

% TEST test_bug1443
% TEST ft_rejectcomponent ft_componentanalysis

load('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat');

cfg=[];
cfg.method='fastica';
cfg.numcomponent = 14; % to make it go fast
cfg.randomseed = 13; % so we get the same output each time
comp = ft_componentanalysis(cfg,data);

assert(ft_senstype(comp,'ctf151'))

cfg = [];
cfg.component = 2; % chosen randomly
rej1 = ft_rejectcomponent(cfg, comp, data);
rej2 = ft_rejectcomponent(cfg, comp);

norm(rej2.grad.tra-rej1.grad.tra)/norm(rej2.grad.tra)
figure;imagesc(rej2.grad.tra - rej1.grad.tra);caxis([-1 1])

load standard_sourcemodel3d10mm
load standard_singleshell
cfg=[];
cfg.grid=sourcemodel;
cfg.vol=vol;

cfg.grad=rej1.grad;
grid1 = ft_prepare_leadfield(cfg, rej1);

cfg.grad=rej2.grad;
grid2 = ft_prepare_leadfield(cfg, rej2);

assert(~isequalwithequalnans(grid1.leadfield,grid2.leadfield))
norm(grid1.leadfield{grid1.inside(1)}-grid2.leadfield{grid2.inside(1)})/norm(grid1.leadfield{grid1.inside(1)})



