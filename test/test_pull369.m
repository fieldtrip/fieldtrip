function test_pull369

% WALLTIME 00:10:00
% MEM 4gb

%%

begsample = (1:100)*1200 - 1199;
endsample = (1:100)*1200;
offset    = (1:100)*0;

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/regressconfound/TacStimRegressConfound.ds');
cfg.trl = [begsample(:) endsample(:) offset(:)];

cfg.numclusters = 1;
grad1 = ft_headmovement(cfg)

%%

cfg.numclusters = 2;
grad2 = ft_headmovement(cfg)

assert(size(grad2.tra,2) == 2*size(grad1.tra,2))

%%

cfg.numclusters = 10;
grad10 = ft_headmovement(cfg)

assert(size(grad10.tra,2) == 10*size(grad1.tra,2))
