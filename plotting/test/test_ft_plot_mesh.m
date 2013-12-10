function test_ft_plot_mesh

% MEM 1gb
% WALLTIME 0:03:00

% TEST test_ft_plot_mesh
% TEST ft_plot_mesh

bnd.pnt = randn(20,3);
bnd.tri = delaunay(bnd.pnt(:,1:2));

figure
ft_plot_mesh(bnd);
