function test_bug1443

load('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat');

cfg=[];
cfg.method='fastica';
cfg.numcomponent = 14; % to make it go fast
cfg.randomseed = 13; % so we get the same output each time
comp = ft_componentanalysis(cfg,data);

cfg = [];
cfg.component = 2; % chosen randomly
rej1 = ft_rejectcomponent(cfg, comp, data);
rej2 = ft_rejectcomponent(cfg, comp);

norm(rej2.grad.tra-rej1.grad.tra)/norm(rej2.grad.tra)
figure;imagesc(rej2.grad.tra - rej1.grad.tra);caxis([-1 1])



