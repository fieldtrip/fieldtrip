function test_bug1149

% MEM 1500mb
% WALLTIME 00:10:00


% generate a unit sphere
[pnt, tri] = icosahedron162;

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


figure,ft_plot_vol(vol,'faceindex','yes')
figure,ft_plot_vol(vol,'faceindex','none')
figure,ft_plot_vol(vol,'facecolor','g')
figure,ft_plot_vol(vol,'facecolor','none')
