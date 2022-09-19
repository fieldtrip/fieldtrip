function inspect_ft_select_range

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_select_range

%%  The following example allows multiple box-like selections to be made
% click in the last selection to trigger the action

close all; clc

x = randn(10,1);
y = randn(10,1);
figure; plot(x, y, '.');

set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'event', 'WindowButtonDownFcn',   'multiple', true, 'xrange', true, 'yrange', true, 'callback', @disp});
set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'event', 'WindowButtonMotionFcn', 'multiple', true, 'xrange', true, 'yrange', true, 'callback', @disp});
set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'event', 'WindowButtonUpFcn',     'multiple', true, 'xrange', true, 'yrange', true, 'callback', @disp});


%% The following example allows a single horizontal selection to be made
close all; clc

x = randn(10,1);
y = randn(10,1);
figure; plot(x, y, '.');

set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'event', 'WindowButtonDownFcn',   'multiple', false, 'xrange', true, 'yrange', false, 'callback', @disp});
set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'event', 'WindowButtonMotionFcn', 'multiple', false, 'xrange', true, 'yrange', false, 'callback', @disp});
set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'event', 'WindowButtonUpFcn',     'multiple', false, 'xrange', true, 'yrange', false, 'callback', @disp});


%% The following example allows multiple horizontal selections to be made
close all; clc

x = randn(10,1);
y = randn(10,1);
figure; plot(x, y, '.');

set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'event', 'WindowButtonDownFcn',   'multiple', true, 'xrange', true, 'yrange', false, 'callback', @disp});
set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'event', 'WindowButtonMotionFcn', 'multiple', true, 'xrange', true, 'yrange', false, 'callback', @disp});
set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'event', 'WindowButtonUpFcn',     'multiple', true, 'xrange', true, 'yrange', false, 'callback', @disp});


%% The following example allows a single vertical selection to be made
close all; clc

x = randn(10,1);
y = randn(10,1);
figure; plot(x, y, '.');

set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'event', 'WindowButtonDownFcn',   'multiple', false, 'xrange', false, 'yrange', true, 'callback', @disp});
set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'event', 'WindowButtonMotionFcn', 'multiple', false, 'xrange', false, 'yrange', true, 'callback', @disp});
set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'event', 'WindowButtonUpFcn',     'multiple', false, 'xrange', false, 'yrange', true, 'callback', @disp});


%%  The following example allows a multiple vertical selections to be made
close all; clc

x = randn(10,1);
y = randn(10,1);
figure; plot(x, y, '.');

set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'event', 'WindowButtonDownFcn',   'multiple', true, 'xrange', false, 'yrange', true, 'callback', @disp});
set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'event', 'WindowButtonMotionFcn', 'multiple', true, 'xrange', false, 'yrange', true, 'callback', @disp});
set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'event', 'WindowButtonUpFcn',     'multiple', true, 'xrange', false, 'yrange', true, 'callback', @disp});


%% The following example allows a single point to be selected
close all; clc

x = 1:10;
y = 1:10;
figure; plot(x, y, '.');

set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'event', 'WindowButtonDownFcn',   'multiple', false, 'xrange', false, 'yrange', false, 'callback', @disp});
set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'event', 'WindowButtonMotionFcn', 'multiple', false, 'xrange', false, 'yrange', false, 'callback', @disp});
set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'event', 'WindowButtonUpFcn',     'multiple', false, 'xrange', false, 'yrange', false, 'callback', @disp});
