function test_bug2560

% WALLTIME 00:15:00
% MEM 300mb

% TEST ft_sourceplot
% TEST ft_plot_slice

% this function is meant to test the correctness of the 'ortho' method in
% sourceplot. It relies on ft_plot_ortho and ft_plot_slice. there may be a
% slight problem in dealing with the functional data/mask

% investigate the issue by creating some data
[x,y,z] = ndgrid(1:7,1:8,1:9);
inside  = zeros(7,8,9);
inside(2:6,2:7,2:8) = 1;
outside = find(inside==0);
inside  = find(inside==1);

source = [];
source.pos = [x(:) y(:) z(:)]; clear x y z
source.dim = [7 8 9];
source.inside  = inside;
source.outside = outside;
source.avg.pow = zeros(7,8,9);
source.avg.pow(inside) = 1;

cfg = [];
cfg.funparameter = 'avg.pow';
cfg.method       = 'ortho';
cfg.funcolormap  = 'jet';
ft_sourceplot(cfg, source);

% the code above confirms an the hunch I had: the inside mask is shifted
% relative to the functional data
