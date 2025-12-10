function test_pull369

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY
% DATA public

%%

begsample = (1:100)*1200 - 1199;
endsample = (1:100)*1200;
offset    = (1:100)*0;

cfg = [];
cfg.dataset = dccnpath('/project/3031000.02/external/download/example/regressconfound/TacStimRegressConfound.ds');
cfg.trl = [begsample(:) endsample(:) offset(:)];

cfg.numclusters = 1;
data1 = ft_headmovement(cfg);

%%

cfg.numclusters = 2;
data2 = ft_headmovement(cfg);

assert(size(data2.grad.tra,2) == 2*size(data1.grad.tra,2))

%%

cfg.numclusters = 10;
data10 = ft_headmovement(cfg);

assert(size(data10.grad.tra,2) == 10*size(data1.grad.tra,2))
