function test_ft_plot_axes

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_plot_axes

close all

result = {};
figure; ft_plot_axes([], 'unit', 'mm'); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm', 'coordsys', 'ctf'); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm', 'coordsys', 'neuromag'); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm', 'coordsys', 'paxinos'); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm', 'coordsys', 'lps'); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm', 'fontcolor', 'y'); view([1 0.5 0.5]); result{end+1} = getframe;
figure; ft_plot_axes([], 'unit', 'mm', 'fontsize', 20); view([1 0.5 0.5]); result{end+1} = getframe;

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}));
  end
end
