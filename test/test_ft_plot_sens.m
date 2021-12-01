function test_ft_plot_sens

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_plot_sens

% make a unit sphere for the head
[headshape.pos, headshape.tri] = mesh_sphere(300);
headshape.unit = 'dm'; % it is 1 dm, 10 cm, or 100 mm

%%

pnt = mesh_sphere(162); % refined icosahedron
pnt(pnt(:,3)<0,:) = []; % only in the upper hemisphere

n = size(pnt,1);
for i=1:n
  pnt(i,:) = pnt(i,:)/norm(pnt(i,:));
  grad.coilpos(i  ,:) = 1.2 * pnt(i,:);
  grad.coilpos(i+n,:) = 1.5 * pnt(i,:);
  grad.coilori(i  ,:) = pnt(i,:);
  grad.coilori(i+n,:) = pnt(i,:);
  grad.tra(i,i)   =  1;
  grad.tra(i,i+n) = -1;
  grad.chanpos(i,:) = 1.02 * pnt(i,:);
  grad.label{i} = num2str(i);
end

% this does not need to be tested exhaustively
grad.coordsys = 'ctf';
axes = true;

figure; ft_plot_sens(grad);
figure; ft_plot_sens(grad, 'axes', axes, 'coil', false);
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true);

figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'coilshape', 'point');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'coilshape', 'circle');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'coilshape', 'square');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'coilshape', 'sphere');

figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', true, 'coilshape', 'point');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', true, 'coilshape', 'circle');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', true, 'coilshape', 'square');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', true, 'coilshape', 'sphere');

figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', false, 'coilshape', 'point');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', false, 'coilshape', 'circle');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', false, 'coilshape', 'square');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', false, 'coilshape', 'sphere');

figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', false, 'orientation', true, 'coilshape', 'point');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', false, 'orientation', true, 'coilshape', 'circle');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', false, 'orientation', true, 'coilshape', 'square');
figure; ft_plot_sens(grad, 'axes', axes, 'coil', true, 'individual', false, 'orientation', true, 'coilshape', 'sphere');

%%

elec.pnt = randn(32,3);
elec.pnt(:,3) = abs(elec.pnt(:,3)); % only in the upper hemisphere
for i=1:32
  elec.pnt(i,:) = elec.pnt(i,:) / norm(elec.pnt(i,:));
  elec.label{i} = num2str(i);
end

% this does not need to be tested exhaustively
elec.coordsys = 'ctf';
axes = true;

figure; ft_plot_sens(elec);
figure; ft_plot_sens(elec, 'axes', axes, 'elec', false);
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true);

figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'elecshape', 'point');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'elecshape', 'circle');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'elecshape', 'square');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'elecshape', 'sphere');

figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', true, 'elecshape', 'point');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', true, 'elecshape', 'circle');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', true, 'elecshape', 'square');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', true, 'elecshape', 'sphere');

figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'elecshape', 'point');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'elecshape', 'circle');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'elecshape', 'square');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'elecshape', 'sphere');

figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'orientation', true, 'elecshape', 'point');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'orientation', true, 'elecshape', 'circle');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'orientation', true, 'elecshape', 'square');
figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'orientation', true, 'elecshape', 'sphere');

% these need the Mathworks computer vision toolbox and do not need to be tested regularly
if false
  figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'orientation', true, 'headshape', headshape, 'elecshape', 'point');
  figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'orientation', true, 'headshape', headshape, 'elecshape', 'circle');
  figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'orientation', true, 'headshape', headshape, 'elecshape', 'square');
  figure; ft_plot_sens(elec, 'axes', axes, 'elec', true, 'individual', false, 'orientation', true, 'headshape', headshape, 'elecshape', 'sphere');
end

%%

opto.optopos = randn(32,3);
opto.optopos(:,3) = abs(opto.optopos(:,3)); % only in the upper hemisphere
for i=1:32
  opto.optopos(i,:) = opto.optopos(i,:) / norm(opto.optopos(i,:));
  if mod(i,2)
    % odd numbers
    opto.label{i} = sprintf('S%02d', ceil(i/2));
  else
    % even numbers
    opto.label{i} = sprintf('D%02d', ceil(i/2));
  end
end

% this does not need to be tested exhaustively
opto.coordsys = 'ctf';
axes = true;

figure; ft_plot_sens(opto);
figure; ft_plot_sens(opto, 'axes', axes, 'opto', false);
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true);

figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'optoshape', 'point');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'optoshape', 'circle');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'optoshape', 'square');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'optoshape', 'sphere');

figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', true, 'optoshape', 'point');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', true, 'optoshape', 'circle');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', true, 'optoshape', 'square');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', true, 'optoshape', 'sphere');

figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'optoshape', 'point');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'optoshape', 'circle');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'optoshape', 'square');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'optoshape', 'sphere');

figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'orientation', true, 'optoshape', 'point');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'orientation', true, 'optoshape', 'circle');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'orientation', true, 'optoshape', 'square');
figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'orientation', true, 'optoshape', 'sphere');

% these need the Mathworks computer vision toolbox and do not need to be tested regularly
if false
  figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'orientation', true, 'headshape', headshape, 'optoshape', 'point');
  figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'orientation', true, 'headshape', headshape, 'optoshape', 'circle');
  figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'orientation', true, 'headshape', headshape, 'optoshape', 'square');
  figure; ft_plot_sens(opto, 'axes', axes, 'opto', true, 'individual', false, 'orientation', true, 'headshape', headshape, 'optoshape', 'sphere');
end
