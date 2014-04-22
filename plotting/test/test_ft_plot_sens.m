function test_ft_plot_sens

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_sens
% TEST ft_plot_sens

elec.pnt = randn(10,3);
for i=1:10
  elec.label{i} = num2str(i);
end

figure
ft_plot_sens(elec);
