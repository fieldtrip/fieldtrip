function test_ft_plot_sens

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_plot_sens

elec.pnt = randn(32,3);
for i=1:32
  elec.pnt(i,:) = elec.pnt(i,:) / norm(elec.pnt(i,:));
  elec.label{i} = num2str(i);
end

figure; ft_plot_sens(elec);
figure; ft_plot_sens(elec, 'elec', false);
figure; ft_plot_sens(elec, 'elec', true);

figure; ft_plot_sens(elec, 'elec', true, 'elecshape', 'point');
figure; ft_plot_sens(elec, 'elec', true, 'elecshape', 'circle');
figure; ft_plot_sens(elec, 'elec', true, 'elecshape', 'square');
figure; ft_plot_sens(elec, 'elec', true, 'elecshape', 'sphere');

figure; ft_plot_sens(elec, 'elec', true, 'individual', true, 'elecshape', 'point');
figure; ft_plot_sens(elec, 'elec', true, 'individual', true, 'elecshape', 'circle');
figure; ft_plot_sens(elec, 'elec', true, 'individual', true, 'elecshape', 'square');
figure; ft_plot_sens(elec, 'elec', true, 'individual', true, 'elecshape', 'sphere');

figure; ft_plot_sens(elec, 'elec', true, 'individual', false, 'elecshape', 'point');
figure; ft_plot_sens(elec, 'elec', true, 'individual', false, 'elecshape', 'circle');
figure; ft_plot_sens(elec, 'elec', true, 'individual', false, 'elecshape', 'square');
figure; ft_plot_sens(elec, 'elec', true, 'individual', false, 'elecshape', 'sphere');

figure; ft_plot_sens(elec, 'elec', true, 'individual', false, 'orientation', true, 'elecshape', 'point');
figure; ft_plot_sens(elec, 'elec', true, 'individual', false, 'orientation', true, 'elecshape', 'circle');
figure; ft_plot_sens(elec, 'elec', true, 'individual', false, 'orientation', true, 'elecshape', 'square');
figure; ft_plot_sens(elec, 'elec', true, 'individual', false, 'orientation', true, 'elecshape', 'sphere');

