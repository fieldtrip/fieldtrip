function test_ft_plot_matrix

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_matrix
% TEST ft_plot_matrix

dat = rand(30,50);
% make first and last couple of columns NaN
dat(:,[1 2 3 48 49 50]) = NaN;

% basic
figure
ft_plot_matrix(dat,'clim', [0 1]);


% hightlightstyle opacity
figure
ft_plot_matrix(dat,'clim',[0 1],'highlightstyle','opacity','highlight',dat>0.5)


% hightlightstyle saturation
figure
ft_plot_matrix(dat,'clim',[0 1],'highlightstyle','saturation','highlight',dat>0.5)
