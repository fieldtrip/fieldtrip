function test_ft_plot_sens

% MEM 1gb
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

% plot electrodes as discs
seg.dim = [50 50 50];
seg.transform = eye(4);
seg.unit = 'mm';
seg.brain = ones(50, 50, 50);
cfg = [];
cfg.tissue = 'brain';
cfg.numvertices = 100;
cfg.spmversion = 'spm12';
mesh = ft_prepare_mesh(cfg, seg);

figure; ft_plot_sens(elec, 'elecshape', 'disc', 'headshape', mesh);
