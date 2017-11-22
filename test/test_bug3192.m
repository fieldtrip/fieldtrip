function test_bug3192

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_plot_mesh, ft_plot_box, ft_plot_vol, ft_plot_dipole, ft_plot_headshape

% one color for all vertex
cfg = [];
cfg.grid.xgrid  = -20:5:20;
cfg.grid.ygrid  = -20:5:20;
cfg.grid.zgrid  = -20:5:20;
grid  = ft_prepare_sourcemodel(cfg);
figure, ft_plot_mesh(grid, 'vertexcolor', 'blue', 'facecolor', 'brain', 'edgecolor', 'skull')

% different colors for each vertex
c = [];
temp = {'r','b'};
for iPos=1:length(grid.pos), c = [c temp{mod(iPos,2)+1}]; end
figure, ft_plot_mesh(grid, 'vertexcolor', c)

figure, ft_plot_box([-1 1 2 3], 'facecolor', 'brain')

% ft_plot_sens, _vol, _headshape, and _dipole just forward to ft_plot_mesh
elecs = [];
elecs.elecpos = [23 42 -31; 69 52 1; 61 67 26; 52 65 45];
elecs.label = {'1' '2' '3' '4'};
figure, ft_plot_sens(elecs, 'edgecolor', 'red')

load(dccnpath('/home/common/matlab/fieldtrip/template/headmodel/standard_singleshell.mat'))
figure, ft_plot_vol(vol, 'edgecolor', 'blue', 'vertexcolor', 'red', 'facecolor', 'brain')

vol.pos = [23 42 -31; 69 52 1; 61 67 26; 52 65 45];
figure, ft_plot_headshape(vol, 'edgecolor', 'blue', 'vertexcolor', 'red', 'facecolor', 'brain')

figure, ft_plot_dipole([1 2 3], [1 2 3], 'color', 'brain')
