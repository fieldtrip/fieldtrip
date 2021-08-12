function test_ft_plot_patch

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_plot_patch

hdat = [1:10 10:-1:1];
vdat = rand(1,10);
vdat = [vdat vdat(end:-1:1)+1];
ft_plot_patch(hdat, vdat)
