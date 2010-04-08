function gui_l (x, y)
%
% This file is needed by FASTICAG

% The load dialog for loading new data
% and new initial guess.

% @(#)$Id: gui_l.m,v 1.4 2004/07/27 13:09:26 jarmo Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables

% Handle to the window
global hf_FastICA_Load;

% Handles to some of the controls in window
global he_FastICA_file;

% What is the load type of load dialog
global g_FastICA_loadType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration options
FIGURENAME = 'FastICA: Load';
FIGURETAG = 'f_FastICALoad';
FIGURESIZE = [x y 450 150];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check to see if this figure is already open - it should not!
% Can't have more than one copy - otherwise the global
% variables and handles can get mixed up.
if ~isempty(findobj('Tag',FIGURETAG))
  error('Error: load dialog already open!');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize some of the controls' values

% What are we loading - who is calling?
caller = get(gcf, 'CurrentObject');

switch get(caller, 'Tag')
 case 'b_LoadData'                             % Do we load new data...
  loadString = 'Load data from variable in Matlab.';
  g_FastICA_loadType = 'data';
  FIGURENAME = 'FastICA: Load data';

 case 'b_LoadGuess'                            % ... or new initial guess?
  loadString = 'Load initial guess for mixing matrix A from variable in Matlab.';
  g_FastICA_loadType = 'guess';
  FIGURENAME = 'FastICA: Load initial guess';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the figure
a = figure('Color',[0.8 0.8 0.8], ...
	   'PaperType','a4letter', ...
	   'Name', FIGURENAME, ...
	   'NumberTitle', 'off', ...
	   'Tag', FIGURETAG, ...
	   'Position', FIGURESIZE, ...
	   'MenuBar', 'none');
set (a, 'Resize', 'off');

hf_FastICA_Load = a;

set(hf_FastICA_Load, 'HandleVisibility', 'callback');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From here on it get's ugly as I have not had time to clean it up


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the frames
pos_l=2;
pos_w=FIGURESIZE(3)-4;
pos_h=FIGURESIZE(4)-4;
pos_t=2;
h_f_load_background = uicontrol('Parent',a, ...
  'BackgroundColor',[0.701961 0.701961 0.701961], ...
  'Position',[pos_l pos_t pos_w pos_h], ...
  'Style','frame', ...
  'Tag','f_load_background');

pos_w=120;
pos_l=FIGURESIZE(3)-(pos_w+2+2);
pos_h=FIGURESIZE(4)-2*4;
pos_t=4;
h_f_load_side = uicontrol('Parent',a, ...
  'BackgroundColor',[0.701961 0.701961 0.701961], ...
  'Position',[pos_l pos_t pos_w pos_h], ...
  'Style','frame', ...
  'Tag','f_load_side');

pos_l=4;
pos_w=FIGURESIZE(3)-8-pos_w-2;
pos_h=FIGURESIZE(4)-8;
pos_t=4;
h_f_load = uicontrol('Parent',a, ...
  'BackgroundColor',[0.701961 0.701961 0.701961], ...
  'Position',[pos_l pos_t pos_w pos_h], ...
  'Style','frame', ...
  'Tag','f_load');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls in f_load
bgc = get(h_f_load, 'BackgroundColor');

pos_w=230;

pos_frame=get(h_f_load, 'Position');
pos_h = 40;
pos_t = pos_frame(2) + pos_frame(4) - pos_h - 6;
pos_l = pos_frame(1) + 6;

b = uicontrol('Parent',a, ...
  'BackgroundColor',bgc, ...
  'HorizontalAlignment','left', ...
  'Position',[pos_l pos_t pos_w pos_h], ...
  'String',loadString, ...
  'Style','text', ...
  'Tag','t_93');

pos_h = 20;
pos_t = pos_t - pos_h - 10;
pos_l = pos_frame(1) + 6;

b = uicontrol('Parent',a, ...
  'BackgroundColor',bgc, ...
  'HorizontalAlignment','left', ...
  'Position',[pos_l pos_t pos_w pos_h], ...
  'String','Name of the variable:', ...
  'Style','text', ...
  'Tag','t_92');

pos_w = 200;
pos_l = pos_l + 30;
pos_t = pos_t - pos_h;
he_FastICA_file = uicontrol('Parent',a, ...
  'BackgroundColor',[1 1 1], ...
  'HorizontalAlignment','left', ...
  'Position',[pos_l pos_t pos_w pos_h], ...
  'String','', ...
  'Style','edit', ...
  'Tag','e_file');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls in f_load_side
pos_vspace = 6;
pos_hspace = 10;
pos_frame = get(h_f_load_side, 'Position');
pos_w = 100;
pos_h = 30;
pos_l = pos_frame(1) + pos_hspace;
pos_t = pos_frame(2) + pos_frame(4) - pos_h - pos_vspace;
b = uicontrol('Parent',a, ...
  'BackgroundColor',[0.701961 0.701961 0.701961], ...
  'Callback','gui_lc Load', ...
  'Position',[pos_l pos_t pos_w pos_h], ...
  'String','Load', ...
  'Tag','b_lLoad');

pos_t=pos_t-pos_h-pos_vspace;
b = uicontrol('Parent',a, ...
  'BackgroundColor',[0.701961 0.701961 0.701961], ...
  'Callback','gui_lc Cancel', ...
  'Position',[pos_l pos_t pos_w pos_h], ...
  'String','Cancel', ...
  'Tag','b_lCancel');

pos_t = pos_frame(2) + pos_vspace;
b = uicontrol('Parent',a, ...
  'BackgroundColor',[0.701961 0.701961 0.701961], ...
  'Callback','gui_lc Help', ...
  'Position',[pos_l pos_t pos_w pos_h], ...
  'String','Help', ...
  'Tag','b_lHelp');

