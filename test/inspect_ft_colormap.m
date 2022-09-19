function inspect_ft_colormap

% see https://github.com/fieldtrip/fieldtrip/pull/1430

close all

dat = linspace(0, 1, 64)';
dat = flipud(dat);

%%

list = {
  'parula'
  'jet'
  'hsv'
  'hot'
  'cool'
  'spring'
  'summer'
  'autumn'
  'winter'
  'gray'
  'bone'
  'copper'
  'pink'
  'lines'
  'colorcube'
  'prism'
  'flag'
  %   'white'
  };

figure
set(gcf, 'Name', 'MATLAB')

n = ceil(sqrt(numel(list)));
m = ceil(numel(list)/n);

for i=1:numel(list)
  subplot(n, m, i);
  imagesc(dat);
  axis off
  disp(list{i})
  % ft_colormap(list{i}); % apply it to the whole figure (which is not interesting)
  ft_colormap(gca, list{i}, 4); % with a very small N
  % ft_colormap(gca, list{i}); % with the default N
  title(list{i})
end

%%

list = {
  'cividis'
  'inferno'
  'magma'
  'plasma'
  'tab10'
  'tab20'
  'tab20b'
  'tab20c'
  'twilight'
  'viridis'
  };

figure
set(gcf, 'Name', 'MATPLOTLIB')

n = ceil(sqrt(numel(list)));
m = ceil(numel(list)/n);

for i=1:numel(list)
  subplot(n, m, i);
  imagesc(dat);
  axis off
  disp(list{i})
  % ft_colormap(list{i}); % apply it to the whole figure (which is not interesting)
  % ft_colormap(gca, list{i}, 4); % with a very small N
  ft_colormap(gca, list{i}); % with the default N
  title(list{i})
end

%%

list = {
  'BrBG'
  'PiYG'
  'PRGn'
  'PuOr'
  'RdBu'
  'RdGy'
  'RdYlBu'
  'RdYlGn'
  'Spectral'
  'Accent'
  'Dark2'
  'Paired'
  'Pastel1'
  'Pastel2'
  'Set1'
  'Set2'
  'Set3'
  'Blues'
  'BuGn'
  'BuPu'
  'GnBu'
  'Greens'
  'Greys'
  'OrRd'
  'Oranges'
  'PuBu'
  'PuBuGn'
  'PuRd'
  'Purples'
  'RdPu'
  'Reds'
  'YlGn'
  'YlGnBu'
  'YlOrBr'
  'YlOrRd'
  };

figure
set(gcf, 'Name', 'BREWERMAP')

n = ceil(sqrt(numel(list)));
m = ceil(numel(list)/n);

for i=1:numel(list)
  subplot(n, m, i);
  imagesc(dat);
  axis off
  disp(list{i})
  % ft_colormap(list{i}); % apply it to the whole figure (which is not interesting)
  % ft_colormap(gca, list{i}, 4); % with a very small N
  ft_colormap(gca, list{i}); % with the default N
  title(list{i})
end


figure
set(gcf, 'Name', 'BREWERMAP REVERSE')

n = ceil(sqrt(numel(list)));
m = ceil(numel(list)/n);

for i=1:numel(list)
  subplot(n, m, i);
  imagesc(dat);
  axis off
  list{i} = ['*' list{i}];
  disp(list{i})
  % ft_colormap(list{i}); % apply it to the whole figure (which is not interesting)
  % ft_colormap(gca, list{i}, 4); % with a very small N
  ft_colormap(gca, list{i}); % with the default N
  title(list{i})
end
