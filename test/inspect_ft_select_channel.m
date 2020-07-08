function inspect_ft_select_channel

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_select_channel

ft_debug off

%% Example 1
close all; clc

% create a figure
figure
cfg = [];
cfg.channel = {'chan1', 'chan2', 'chan3', 'chan4', 'chan5'};
cfg.layout  = 'ordered';
cfg.rows = 3;
cfg.columns = 3;
lay = ft_prepare_layout(cfg);
ft_plot_layout(lay)

% add the required guidata
info       = guidata(gcf);
info.x     = lay.pos(:,1);
info.y     = lay.pos(:,2);
info.label = lay.label;
guidata(gcf, info)

% add this function as the callback to make a single selection
set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'callback', @disp})

% or to make multiple selections
set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', false, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', false, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', false, 'callback', @disp, 'event', 'WindowButtonDownFcn'})

%% Example 2, executed from within a subplot
close all; clc

% create a figure
figure
subplot(2,2,1)
cfg = [];
cfg.channel = {'chan1', 'chan2', 'chan3', 'chan4'};
cfg.layout  = 'ordered';
lay = ft_prepare_layout(cfg);
ft_plot_layout(lay)

% add the channel information to guidata under identifier linked to this axis
ident              = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
set(gca,'tag',ident);
info               = guidata(gcf);
info.(ident).x     = lay.pos(:, 1);
info.(ident).y     = lay.pos(:, 2);
info.(ident).label = lay.label;
guidata(gcf, info);

% add this function as the callback to make a single selection
set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'callback', @disp})

% or to make multiple selections
set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
