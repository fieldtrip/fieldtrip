function test_ft_plot_cloud

% WALLTIME 00:10:00
% MEM 2gb

elec = ft_read_sens('easycap-M10.txt');

pos = elec.chanpos;
val = linspace(-1, 1, numel(elec.label));  % values between -1 and 1

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




