function inspect_ft_select_box

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_select_box

ft_debug off

%% Example 1
close all; clc

x = rand(1, 100);
y = rand(1, 100);

figure
a1 = subplot(1,2,1); plot(-x, -y, '.')
a2 = subplot(1,2,2); plot(+x, +y, '.')

[x, y] = ft_select_box()

