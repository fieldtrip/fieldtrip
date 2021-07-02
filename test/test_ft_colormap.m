function test_ft_colormap

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_colormap

cmap = {'parula', 'jet', 'hsv', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'lines', 'colorcube', 'prism', 'flag', 'cividis', 'inferno', 'magma', 'plasma', 'tab10', 'tab20', 'tab20b', 'tab20c', 'twilight', 'viridis', 'thermal', 'haline', 'solar', 'ice', 'oxy', 'deep', 'dense' 'algae', 'matter', 'turbid', 'speed', 'amp', 'tempo', 'rain', 'delta', 'curl', 'diff', 'tarn', 'phase', 'topo'};

n = 100;

close all
figure
imagesc(1:n)

for i=1:length(cmap)
  % test the colormap
  disp(cmap{i});
  ft_colormap(cmap{i}, n);
  drawnow
  % test the flipped colormap by adding a minus in front of it
  disp(['-' cmap{i}]);
  ft_colormap(['-' cmap{i}], n);
  drawnow
end
