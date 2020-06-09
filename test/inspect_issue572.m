function inspect_issue572

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

close all

data.label = {'test'};
data.time = {1:1000};
data.trial = {sin(1:1000)};
data.fsample = 1000;
data.sampleinfo = [1 1000];

cfg = [];
ft_databrowser(cfg, data);  % plot fake data

figure('Position', [1 1 100 100]);  % open a new figure

% when the mouse cursor passes over the ft_databrowser figure, I get the
% following errors, continuously appearing in the command window until I
% click anywhere within the ft_databrowser figure:

% Error in ft_select_range (line 147)
%   set(findobj(handle,'hittest','on'), 'uicontextmenu',hcmenu);
%
% Error while evaluating Figure WindowButtonMotionFcn
%
% Error using matlab.graphics.Graphics/set
% The UIContextMenu's parent figure is not the same as the parent figure of this object
