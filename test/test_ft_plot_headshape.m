function test_ft_plot_headshape

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_headshape
% TEST ft_plot_headshape

shape.pnt = randn(500,3);
shape.fid.pnt = randn(3,3);
shape.fid.label = {'nas', 'lpa', 'rpa'};

figure
ft_plot_headshape(shape);
