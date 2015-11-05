function test_ft_plot_matrix

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_matrix
% TEST ft_plot_matrix

dat  = rand(30,50);
mask = rand(30,50);
% make first and last couple of columns NaN
dat(:,[1 2 3 48 49 50])  = NaN;
mask(:,[1 2 3 48 49 50]) = 0; % make the mask also a NaN mask

% 1 basic
figure
ft_plot_matrix(dat,'clim', [0 1]);


% 2 hightlightstyle opacity binary mask
figure
ft_plot_matrix(dat,'clim',[0 1],'highlightstyle','opacity','highlight',mask>0.5)


% 3 hightlightstyle saturation binary mask
figure
ft_plot_matrix(dat,'clim',[0 1],'highlightstyle','saturation','highlight',mask>0.5)


% 4 hightlightstyle opacity graded mask
figure
ft_plot_matrix(dat,'clim',[0 1],'highlightstyle','opacity','highlight',mask)


% 5 hightlightstyle saturation graded mask
figure
ft_plot_matrix(dat,'clim',[0 1],'highlightstyle','saturation','highlight',mask)
