function inspect_bug3325

% WALLTIME 00:30:00
% MEM 1gb

% this option makes it interactive
interactive = false;

%%

grad = [];
grad.label = {'1', '2', '3'};
grad.chantype = {'meggrad', 'meggrad', 'megmag'};
grad.unit = 'mm';
grad.coilpos = [
  0 0  0
  40 0  0
  80 0  0
  0 0 20
  40 0 20
  ];
grad.coilori = [
  0 0 1
  0 0 1
  0 0 1
  0 0 1
  0 0 1
  ];
grad.tra = [
  1 0 0 -1  0
  0 1 0  0 -1
  0 0 1  0  0
  ];
grad.chanpos = grad.coilpos([1 2 5],:);
grad.chanori = grad.coilori([1 2 5],:);

%%

arg = [];
arg(1).name = 'coil';
arg(1).value = {true, false};
arg(2).name = 'chantype';
arg(2).value = {[], 'meggrad', 'megmag'};
arg(3).name = 'coilshape';
arg(3).value = {'point', 'circle', 'square', 'sphere'};
arg(4).name = 'marker';
arg(4).value = {'.', '*'}; % transpose('.*+')
arg(5).name = 'style';
arg(5).value = {'r', 'r^' 'bv'};  % ['r'; 'g'; 'b'], ['r.'; 'b^'; 'gv']
arg(6).name = 'fontcolor';
arg(6).value = {'k', 'skin'};
arg(7).name = 'facecolor';
arg(7).value = {'g', 'none'};
arg(8).name = 'edgecolor';
arg(8).value = {'b', 'red'};


opt = constructalloptions(arg);

%%

close all

for i=1:size(opt,1)
  figure
  ft_plot_sens(grad, opt{i,:});
  view(20,20)
  % drawnow
  fprintf('-------------------------i = %d -------------------------\n', i);
  disp(opt(i,:))
  if interactive
    disp('press a key')
    pause
  end
  close
end
