function test_bug1149

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY


% generate a unit sphere
[pnt, tri] = mesh_sphere(162);

% create the sphere
bnd.pnt = pnt * 90;
bnd.tri = tri;

vol = [];
vol.bnd = bnd;
vol.type = 'dipoli';

figure,ft_plot_mesh(bnd,'faceindex','yes')
figure,ft_plot_mesh(bnd,'faceindex','none')
figure,ft_plot_mesh(bnd,'facecolor','g')
figure,ft_plot_mesh(bnd,'facecolor','none')


figure,ft_plot_headmodel(vol,'faceindex','yes')
figure,ft_plot_headmodel(vol,'faceindex','none')
figure,ft_plot_headmodel(vol,'facecolor','g')
figure,ft_plot_headmodel(vol,'facecolor','none')
