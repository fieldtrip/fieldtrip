function test_ft_plot_crosshair

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_plot_crosshair

close all

result = {};
figure; ft_plot_axes([], 'unit', 'mm'); ft_plot_crosshair([100 100 100], 'color', 'm'); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm'); ft_plot_crosshair([100 100 100], 'color', 'y'); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm'); ft_plot_crosshair([100 100 100], 'color', [0 0 0]); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm'); ft_plot_crosshair([100 100 100], 'color', 'powderblue'); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm'); ft_plot_crosshair([100 100 100], 'color', 'limegreen'); view([1 0.5 0.5]); result{end+1} = getframe;

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}));
  end
end
