function test_ft_plot_montage

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_plot_montage

close all

n = 101;
img = rand(n,n,n);
img(1,:,:) = 2;
img(:,1,:) = 2;
img(:,:,1) = 2;
img(n,:,:) = 2;
img(:,n,:) = 2;
img(:,:,n) = 2;

transform = eye(4);

result = {};
figure; ft_plot_montage(img, 'transform', transform, 'orientation', [0 0 1]); result{end+1} = getframe;
figure; ft_plot_montage(img, 'transform', transform, 'orientation', [0 1 0]); result{end+1} = getframe;
figure; ft_plot_montage(img, 'transform', transform, 'orientation', [1 0 0]); result{end+1} = getframe;
figure; ft_plot_montage(img, 'transform', transform, 'orientation', [1 0 0], 'nslice', 16); result{end+1} = getframe;
figure; ft_plot_montage(img, 'transform', transform, 'orientation', [1 0 0], 'nslice', 25); result{end+1} = getframe;

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}));
  end
end
