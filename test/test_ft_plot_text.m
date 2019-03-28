function test_ft_plot_text

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_plot_text
% TEST ft_plot_text

%%

figure
ft_plot_text(0.5, 0.5, 'test');

%%
% should be in center, i.e. at [1 1] relative to MATLAB axes

figure
ft_plot_text(0, 0 , '+',    'hpos', 1, 'vpos', 1, 'width', 1, 'height', 1, 'hlim', [-1 1], 'vlim', [-1 1])
axis([0 2 0 2])

%%
% should be in center of the local axes

figure
x = -1:0.1:1;
y = -1:0.1:1;
ft_plot_vector(x, y,        'hpos', 1, 'vpos', 1, 'width', 0.2, 'height', 0.2, 'box', true)
ft_plot_text(0, 0 , '+',    'hpos', 1, 'vpos', 1, 'width', 0.2, 'height', 0.2, 'hlim', [-1 1], 'vlim', [-1 1])
axis([0 2 0 2])

%%
% should be in lower left corner of the local axes

figure
x = 0:0.1:1;
y = 0:0.1:1;
ft_plot_vector(x, y,        'hpos', 1, 'vpos', 1, 'width', 0.2, 'height', 0.2, 'box', true)
ft_plot_text(0, 0 , '+',    'hpos', 1, 'vpos', 1, 'width', 0.2, 'height', 0.2, 'hlim', [0 1], 'vlim', [0 1])
axis([0 2 0 2])
