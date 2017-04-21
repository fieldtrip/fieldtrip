function test_ft_plot_cloud

% WALLTIME 00:10:00
% MEM 2gb

elec = ft_read_sens('easycap-M10.txt');

pos = elec.chanpos;
val = linspace(-1, 1, numel(elec.label));  % values between -1 and 1

% Create test meshes
seg.dim = [50 50 50];
seg.transform = eye(4);
seg.unit = 'mm';
seg.brain = ones(50, 50, 50);
cfg = [];
cfg.tissue = 'brain';
cfg.numvertices = 100;
mesh1 = ft_prepare_mesh(cfg, seg);

seg.brain = zeros(50, 50, 50);
seg.brain(11:40, 11:40, 11:40) = ones(30, 30, 30);
mesh2 = ft_prepare_mesh(cfg, seg);

%%
% the following three figures should be identical, since the units should be correctly detected

% in mm
figure; ft_plot_cloud(pos, val);
% in cm
figure; ft_plot_cloud(pos/10, val);
% in m
figure; ft_plot_cloud(pos/1000, val);

%%
% these should be the same, except for the axes

figure; ft_plot_cloud(pos, val, 'unit', 'mm'); axis on; grid on; xlabel('mm')
figure; ft_plot_cloud(pos, val, 'unit', 'cm'); axis on; grid on; xlabel('cm')
figure; ft_plot_cloud(pos, val, 'unit', 'm');  axis on; grid on; xlabel('m')

%%
% 1 and 2 should be the same, 3 is different

figure; ft_plot_cloud(pos, val, 'colormap', 'default');
figure; ft_plot_cloud(pos, val, 'colormap', 'parula');
figure; ft_plot_cloud(pos, val, 'colormap', 'jet');


%%
% these should all work quite the same

figure; ft_plot_topo(pos(:,1), pos(:,2), val); colorbar
figure; ft_plot_topo3d(pos, val); colorbar
figure; ft_plot_cloud(pos, val); colorbar

%%
% go over the other options

figure; ft_plot_cloud(pos, val, 'ptdensity', 5);
figure; ft_plot_cloud(pos, val, 'ptdensity', 50);
figure; ft_plot_cloud(pos, val, 'colorgrad', 1);
figure; ft_plot_cloud(pos, val, 'scalerad', 'no');

%% 
% mesh options
figure; ft_plot_cloud(pos, val, 'mesh', mesh1);
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'facealpha', [.5 1]); % should be able to see one mesh inside the other

% these two should be the same
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'facealpha', [.5 1], 'facecolor', {'r', 'b'});
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'facealpha', [.5 1], 'facecolor', {[1 0 0], [0 0 1]});

% these two should be the same
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'facealpha', 0, 'edgealpha', [.5 1], 'edgecolor', {'r', 'b'});
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'facealpha', 0, 'edgealpha', [.5 1], 'edgecolor', {[1 0 0], [0 0 1]});

%%
% slice options
% default 2d
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d');

% 2d and 3d should cut the same slices
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '3d', 'nslices', 3);
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d', 'nslices', 3);

% increasing minspace should spread the slices out
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '3d', 'nslices', 3, 'minspace', 20);
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d', 'nslices', 3, 'minspace', 20);

% test 'slicepos'
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '3d', 'slicepos', 33);
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d', 'slicepos', 33);

% test 'slicetype' = 'surf'
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d', 'slicepos', 33, 'slicetype', 'surf');

% test slicing through multiple meshes 
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'slice', '2d', 'slicepos', 33) 

% test slicing in other orientations
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '3d', 'ori', 'x', 'slicepos', 33);
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d', 'ori', 'x', 'slicepos', 33);
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '3d', 'ori', 'z', 'slicepos', 30);
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d', 'ori', 'z', 'slicepos', 30);

% test intersect line settings
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d', 'slicepos', 33, 'intersectcolor', 'r');
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d', 'slicepos', 33, 'intersectlinestyle', '--');
figure; ft_plot_cloud(pos, val, 'mesh', mesh1, 'slice', '2d', 'slicepos', 33, 'intersectlinewidth', 10);
% should allow for line setting inputs for each individual mesh, but default 
% to the first one if there is only 1
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'slice', '2d', 'slicepos', 33, 'intersectcolor', 'r');
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'slice', '2d', 'slicepos', 33, 'intersectcolor', {[0 1 0], [1 0 0]});
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'slice', '2d', 'slicepos', 33, 'intersectlinestyle', '--');
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'slice', '2d', 'slicepos', 33, 'intersectlinestyle', {'--', '-'});
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'slice', '2d', 'slicepos', 33, 'intersectlinewidth', 10);
figure; ft_plot_cloud(pos, val, 'mesh', {mesh1, mesh2}, 'slice', '2d', 'slicepos', 33, 'intersectlinewidth', [10 5]);


