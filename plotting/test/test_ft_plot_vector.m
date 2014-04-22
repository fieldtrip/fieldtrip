function test_ft_plot_vector

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_vector
% TEST ft_plot_vector

mask = zeros(1, 100)
%boundary conditions
mask(1:10) = 1;
mask(90:100) = 1;
mask(40:60) = 1;

figure
subplot(2,2,1)
ft_plot_vector(randn(1,100), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0);
subplot(2,2,2)
ft_plot_vector(randn(1,100), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0, 'highlight', mask);
subplot(2,2,3)
ft_plot_vector(randn(1,100), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0, 'highlight', mask, 'highlightstyle', 'saturation');
subplot(2,2,4)
ft_plot_vector(randn(1,100), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0, 'highlight', mask, 'highlightstyle', 'thickness');
