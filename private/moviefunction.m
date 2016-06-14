function moviefunction(cfg, varargin)

% we need cfg.plotfun to plot the data
% data needs to be 3D, N x time x freq (last can be singleton)
%   N needs to correspond to number of vertices (channels, gridpoints, etc)

% new UI artwork
%
% [main window] -------------------------------------\
% | [uipanel: plot]            | [uipanel: colormap] | 
% |  <80%> ^70%v               |   <20%>             |
% |  [axes: mainaxes]          |  [axes: coloraxes]  |
% |   [axes: subaxes]          |                     |  
% |                            |                     |
% |                            |                     |
% |                            |                     | 
% |                            |                     |
% |                            |                     |
% |                            |                     |  
% |                            |                     |
% |                            |                     |
% |                            |                     |
% |                            |                     |
% |                            |[Checkbox: automatic]|
% |                            |[Checkbox: symmetric]|
% |--[uipanel: control] <100%>-----------------------|   
% | [Button:Go] [Checkbox:Rec]                       |
% | [Slider: frequency]                              |
% | [Slider: time]                                   |
% | [TextField: speed]                               |
% \--------------------------------------------------/
%
% Right click on mainaxes: 
%     rotation options (+automatic), zooming options, open subplot (new subaxes added)
%     start, stop, record, time, freq, colorbar
% Left click on subaxes:  dragging it around
% Right click on subaxes: close
%
% Consider following neat extras
%   new figure for zooming and and rotation timing options during playback
%   add an intro a la Jan-Mathijs for source, similar for normal topos
%   make panels foldable (like on mobile devices)

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim',         'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'renamed',    {'zparam',       'funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'parameter',    'funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'mask',         'maskparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'framespersec', 'framerate'});
cfg = ft_checkconfig(cfg, 'deprecated', {'xparam'});
cfg = ft_checkconfig(cfg, 'forbidden',  {'yparam'});

% set default cfg settings
cfg.xlim            = ft_getopt(cfg, 'xlim',   'maxmin');
cfg.ylim            = ft_getopt(cfg, 'ylim',   'maxmin');
cfg.zlim            = ft_getopt(cfg, 'zlim',   'maxmin');
cfg.xparam          = ft_getopt(cfg, 'xparam', 'time');
cfg.yparam          = ft_getopt(cfg, 'yparam', []);
cfg.maskparameter   = ft_getopt(cfg, 'maskparameter');
cfg.inputfile       = ft_getopt(cfg, 'inputfile',    []);
cfg.moviefreq       = ft_getopt(cfg, 'moviefreq',    []);
cfg.movietime       = ft_getopt(cfg, 'movietime',    []);
cfg.movierpt        = ft_getopt(cfg, 'movierpt',     1);

for i=1:numel(varargin)
  if ~isempty(cfg.yparam) && ~isfield(varargin{i}, 'freq')
    error('data argument %i does not have a .freq field, while all former had', i);
  elseif isempty(cfg.yparam) && isfield(varargin{i}, 'freq')
    cfg.yparam        = 'freq';
  else
    cfg.yparam        = [];
  end
end

% set guidata flags
opt                 = [];
opt.valx            = 1;
opt.valy            = 1;
opt.valz            = 1;
opt.record          = ~istrue(ft_getopt(cfg, 'interactive', 'yes'));
opt.framesfile      = ft_getopt(cfg, 'framesfile',   []);
opt.samperframe     = ft_getopt(cfg, 'samperframe',  1);
opt.framerate       = ft_getopt(cfg, 'framerate',    5);          
opt.fixedframesfile = ~isempty(opt.framesfile);
opt.xvalues   = varargin{1}.(cfg.xparam); % below consistency is checked
if ~isempty(cfg.yparam)
  opt.yvalues = varargin{i}.(cfg.yparam); % below consistency is checked
else
  opt.yvalues   = [];
end

% check data options and consistency
Ndata = numel(varargin);
for i = 1:Ndata
  opt.issource(i)   = ft_datatype(varargin{i}, 'source');
  opt.isfreq(i)     = ft_datatype(varargin{i}, 'freq');
  opt.istimelock(i) = ft_datatype(varargin{i}, 'timelock');
  if opt.issource(i)+opt.isfreq(i)+opt.istimelock(i) ~= 1
    error('data argument %i cannot be definitely identified as frequency, timelock or source data', i);
  end
end

if Ndata>1,
  if all(opt.issource==1),
    % this is allowed maybe in the future: multiple source input arguments
    error('currently, not more than 1 data argument is allowed, unless one of them is a parcellation, and the other a channel level structure');
  elseif all(opt.isfreq==1),
    % this is allowed maybe in the future: multiple freq input arguments
    error('currently, not more than 1 data argument is allowed, unless one of them is a parcellation, and the other a channel level structure');
  elseif all(opt.istimelock==1),
    % this is allowed maybe in the future: multiple timelock input arguments
    error('currently, not more than 1 data argument is allowed, unless one of them is a parcellation, and the other a channel level structure');
  elseif Ndata==2 && sum(opt.issource==1)
    % this may be allowed, provided the source data is a parcellation and
    % the other argument a parcellated data structure
    if ~ft_datatype(varargin{opt.issource}, 'parcellation')
      error('the source data structure should be a parcellation');
    end
  end
end

% set the funparameter
if any(opt.isfreq)
  opt.funparameter = ft_getopt(cfg, 'funparameter', 'powspctrm'); % use power as default
elseif any(opt.istimelock)
  opt.funparameter = ft_getopt(cfg, 'funparameter', 'avg'); % use power as default
elseif any(opt.issource)
  % FIXME a call to ft_checkdata wouldn't hurt here :-)
  opt.funparameter = ft_getopt(cfg, 'funparameter', 'avg.pow'); % use power as default
end

if any(opt.issource) && (any(opt.isfreq) || any(opt.istimelock))
  opt.anatomy  = varargin{opt.issource};
  opt.ismesh   = isfield(opt.anatomy, 'tri');
  opt.isvolume = isfield(opt.anatomy, 'dim');
  
  varargin     = varargin(~opt.issource);
  opt.dat{1}   = getsubfield(varargin{1}, opt.funparameter);
  funtok = tokenize(opt.funparameter, '.');
  opt.dimord   = getdimord(varargin{1}, funtok{end});
elseif any(opt.issource)
  opt.anatomy  = varargin{opt.issource};
  opt.ismesh   = isfield(opt.anatomy, 'tri');
  opt.isvolume = isfield(opt.anatomy, 'dim');
  
  opt.dat{1}   = getsubfield(varargin{opt.issource}, opt.funparameter);
  funtok = tokenize(opt.funparameter, '.');
  opt.dimord   = getdimord(varargin{1}, funtok{end});
else
  opt.layout   = ft_prepare_layout(cfg); % let ft_prepare_layout do the error handling for now
  
  % ensure that the data in all inputs has the same channels, time-axis, etc.
  tmpcfg = keepfields(cfg, {'frequency', 'latency', 'channel'});
  [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
  % restore the provenance information
  [cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

  [opt.seldat, opt.sellay] = match_str(varargin{1}.label, opt.layout.label);
  opt.chanx = opt.layout.pos(opt.sellay,1);
  opt.chany = opt.layout.pos(opt.sellay,2);
  
  for i = 1:numel(varargin)
    opt.dat{i} = getsubfield(varargin{i}, opt.funparameter);
  end
  opt.dimord = getdimord(varargin{1}, opt.funparameter);
end

dimtok = tokenize(opt.dimord, '_');
if any(strcmp(dimtok, 'rpt') | strcmp(dimtok, 'subj'))
  error('the input data cannot contain trials or subjects, please average first using ft_selectdata');
end

opt.ydim = find(strcmp(cfg.yparam, dimtok));
opt.xdim = find(strcmp(cfg.xparam, dimtok));
opt.zdim = setdiff(1:ndims(opt.dat{1}), [opt.ydim opt.xdim]);
if opt.zdim ~=1
  error('input %i data does not have the correct format of N x time (x freq)', i);
end

% permute the data matrix
for i = 1:numel(opt.dat)
  opt.dat{i} = permute(opt.dat{i}, [opt.zdim(:)' opt.xdim opt.ydim]);
end  
opt.xdim = 2;
if ~isempty(cfg.yparam)
  opt.ydim = 3;
else
  opt.ydim = [];
end

% get the mask
for i = 1:numel(varargin)
  if ~isempty(cfg.maskparameter) && ischar(cfg.maskparameter)
    opt.mask{i} = double(getsubfield(varargin{i}, cfg.maskparameter)~=0);
  else
    opt.mask{i} = ones(size(opt.dat{i}));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw the figure with the GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt = createGUI(opt);

if isempty(cfg.yparam)
%  set(opt.handles.label.yparam, 'Visible', 'off');
%  set(opt.handles.slider.yparam, 'Visible', 'off');
%  opt.yparam = '';
else
  opt.yparam = [upper(cfg.yparam(1)) lower(cfg.yparam(2:end))];
  set(opt.handles.label.yparam, 'String', opt.yparam);
end
opt.xparam = [upper(cfg.xparam(1)) lower(cfg.xparam(2:end))];
set(opt.handles.label.xparam, 'String', opt.xparam);

% set timer
opt.timer = timer;
set(opt.timer, 'timerfcn', {@cb_timer, opt.handles.figure}, 'period', 0.1, 'executionmode', 'fixedSpacing');

opt = plot_geometry(opt);
opt = plot_other(opt);

opt.speed = 1;
guidata(opt.handles.figure, opt);
% init first screen

cb_slider(opt.handles.figure);
%uicontrol(opt.handles.colorbar);
if strcmp(cfg.zlim, 'maxmin') 
  set(opt.handles.checkbox.automatic, 'Value', 1) 
  set(opt.handles.checkbox.symmetric, 'Value', 0)
  cb_colorbar(opt.handles.figure);
elseif strcmp(cfg.zlim, 'maxabs')
  set(opt.handles.checkbox.automatic, 'Value', 1)
  set(opt.handles.checkbox.symmetric, 'Value', 1)
  cb_colorbar(opt.handles.figure);
else
  set(opt.handles.checkbox.automatic, 'Value', 0)
  set(opt.handles.checkbox.symmetric, 'Value', 0)
  cb_colorbar(opt.handles.figure, cfg.zlim);
end

if opt.record
  opt.record = false;
  opt.quit = true;
  guidata(opt.handles.figure, opt);
  cb_recordbutton(opt.handles.button.record);
else 
  opt.quit = false;
  guidata(opt.handles.figure, opt);
end

end

%%
%  **************************************************************
%  ********************* CREATE GUI *****************************
%  **************************************************************
function opt = createGUI(opt)

%% main figure
opt.handles.figure = figure(    ...
  'Units',       'normalized',  ...
  'Name',        'ft_movieplot',...
  'Menu',        'none',        ...
  'Toolbar',     'figure',      ...
  'NumberTitle', 'off',         ...
  'UserData',    [],            ...
  'Tag',         'mainFigure',  ...
  'Visible',     'on',          ...
  'windowbuttondownfcn', @cb_getposition);
  %'WindowButtonUpFcn', @cb_stopDrag, ...

%%  panels 
clf

% visualization panel for the geometrical information
opt.handles.panel.visualization_geometry = uipanel(...
  'tag',      'mainPanels1',      ... % tag according to position
  'parent',   opt.handles.figure, ...
  'units',    'normalized',       ...
  'title',    'Visualization_geometry', ...
  'clipping', 'on',               ...
  'visible',  'on');

% rearrange 
ft_uilayout(opt.handles.figure,     ...
  'tag',             'mainPanels1', ...
  'backgroundcolor', [.8 .8 .8],    ...
  'hpos',            'auto',        ...
  'vpos',            .15,           ...
  'halign',          'left',        ...
  'width',           1,             ...
  'height',          0.1);

% visualization panel for the non-geometrical information
opt.handles.panel.visualization = uipanel(...
  'tag',      'mainPanels2',      ... % tag according to position
  'parent',   opt.handles.figure, ...
  'units',    'normalized',       ...
  'title',    'Visualization',    ...
  'clipping', 'on',               ...
  'visible',  'on');

% rearrange 
ft_uilayout(opt.handles.figure,     ...
  'tag',             'mainPanels2', ...
  'backgroundcolor', [.8 .8 .8],    ...
  'hpos',            'auto',        ...
  'vpos',            .55,           ...
  'halign',          'left',        ...
  'width',           0.5,           ...
  'height',          0.1);

% settings panel (between different views can be switched)
% opt.handles.panel.view = uipanel(...
%   'tag',    'sidePanels', ... % tag according to position
%   'parent', opt.handles.figure,...
%   'units',  'normalized', ...
%   'title',  'View',...
%   'clipping','on',...
%   'visible', 'on');

% settings panel ()
opt.handles.panel.settings = uipanel(...
  'tag',      'sidePanels',       ... % tag according to position
  'parent',   opt.handles.figure, ...
  'units',    'normalized',       ...
  'title',    'Settings',         ...
  'clipping', 'on',               ...
  'visible',  'on');

% rearrange panel
ft_uilayout(opt.handles.figure,    ...
  'tag',             'sidePanels', ...
  'backgroundcolor', [.8 .8 .8],   ...
  'hpos',            'auto',       ...
  'vpos',            0.15,         ...
  'halign',          'right',      ...
  'width',           0.15,         ...
  'height',          0.55);

% control panel
opt.handles.panel.controls = uipanel(...
  'tag',    'lowerPanels', ... % tag according to position
  'parent', opt.handles.figure,...
  'units',  'normalized', ...
  'title',  'Controls',...
  'clipping','on',...
  'visible', 'on');

% rearrange 
ft_uilayout(opt.handles.figure,     ...
  'tag',             'lowerPanels', ...
  'backgroundcolor', [0.8 0.8 0.8], ...
  'hpos',            'auto',        ...
  'vpos',            0,             ...
  'halign',          'right',       ...
  'width',           1,             ...
  'height',          0.15);

ft_uilayout(opt.handles.figure, ...
  'tag',   'mainPanels1',       ...
  'retag', 'sidePanels');

ft_uilayout(opt.handles.figure, ...
  'tag',   'mainPanels2',       ...
  'retag', 'sidePanels');

ft_uilayout(opt.handles.figure,    ...
  'tag',             'sidePanels', ...
  'backgroundcolor', [.8 .8 .8],   ...
  'hpos',            'auto',       ...
  'vpos',            'align',      ...
  'halign',          'left',       ...
  'valign',          'top',        ...
  'height',          .85);


% add axes
% 3 axes for switching to topo viewmode
% opt.handles.axes.A = axes(...
%   'Parent',opt.handles.panel.view,...
%   'Position',[0 0 1 .33], ...
%   'Tag','topoAxes',...
%   'HandleVisibility', 'on', ...
%   'UserData',[]);
% 
% opt.handles.axes.B = axes(...
%   'Parent',opt.handles.panel.view,...
%   'Position',[0 0.33 1 .33], ...
%   'Tag','topoAxes',...
%   'HandleVisibility', 'on', ...
%   'UserData',[]);
% 
% opt.handles.axes.C = axes(...
%   'Parent',opt.handles.panel.view,...
%   'Position',[0 0.66 1 .33], ...
%   'Tag','topoAxes',...
%   'HandleVisibility', 'on', ...
%   'UserData',[]);

% control panel uicontrols

opt.handles.button.play = uicontrol(...
  'parent',opt.handles.panel.controls,...
  'tag', 'controlInput', ...
  'string','Play',... 
  'units', 'normalized', ...
  'style', 'pushbutton', ...
  'callback',@cb_playbutton,...
  'userdata', 'space');
  
opt.handles.button.record = uicontrol(...
  'parent',opt.handles.panel.controls,...
  'tag', 'controlInput', ...
  'string','Record?',... 
  'units', 'normalized', ...
  'style', 'checkbox', ...
  'Callback',@cb_recordbutton,...
  'userdata', 'enter');

ft_uilayout(opt.handles.panel.controls, ...
  'tag', 'controlInput', ...
  'width', 0.15, ...
  'height', 0.25, ...
  'hpos', 'auto', ...
  'vpos', .75, ...
  'halign', 'right', ...
  'valign', 'top');


opt.handles.text.speed = uicontrol(...
  'parent',opt.handles.panel.controls,...
  'tag', 'speedInput', ...
  'string','Speed',... 
  'units', 'normalized', ...
  'style', 'text');

opt.handles.edit.speed = uicontrol(...
  'parent',opt.handles.panel.controls,...
  'tag', 'speedInput', ...
  'string','1.0',... 
  'units', 'normalized', ...
  'Callback',@cb_speed, ...
  'style', 'edit');


ft_uilayout(opt.handles.panel.controls, ...
  'tag', 'speedInput', ...
  'width', 0.15, ...
  'height', 0.25, ...
  'hpos', 'auto', ...
  'vpos', .25, ...
  'halign', 'right', ...
  'valign', 'top');


% opt.handles.button.faster = uicontrol(...
%   'parent',opt.handles.panel.controls,...
%   'tag', 'controlInput', ...
%   'string','+',... 
%   'units', 'normalized', ...
%   'style', 'pushbutton', ...
%   'userdata', '+', ...
%   'Callback',@cb_slider);
% 
% 
% opt.handles.button.slower = uicontrol(...
%   'parent',opt.handles.panel.controls,...
%   'tag', 'controlInput', ...
%   'string','-',... 
%   'units', 'normalized', ...
%   'style', 'pushbutton', ...
%   'userdata', '-', ...
%   'Callback',@cb_slider);

% ft_uilayout(opt.handles.panel.controls, ...
%   'tag', 'controlInput', ...
%   'width', 0.15, ...
%   'height', 0.25, ...
%   'hpos', 'auto', ...
%   'vpos', .75, ...
%   'valign', 'top');
% 
% set(opt.handles.button.slower, 'tag', 'speedButtons');
% set(opt.handles.button.faster, 'tag', 'speedButtons');
% 
% 
% ft_uilayout(opt.handles.panel.controls, ...
%   'tag', 'speedButtons', ...
%   'width', 0.05, ...
%   'height', 0.125, ...
%   'vpos', 'auto');
% 
% ft_uilayout(opt.handles.panel.controls, ...
%   'tag', 'speedButtons', ...
%   'hpos', 'align');
% 
% 
% ft_uilayout(opt.handles.panel.controls, ...
%   'tag', 'speedButtons', ...
%   'retag', 'controlInput');

% speed control

opt.handles.slider.xparam = uicontrol(...
  'Parent',opt.handles.panel.controls,...
  'Units','normalized',...
  'BackgroundColor',[0.9 0.9 0.9],...
  'Callback',@cb_slider,...
  'Position',[0.0 0.75 0.625 0.25],...
  'SliderStep', [1/numel(opt.xvalues) 1/numel(opt.xvalues)],...  
  'String',{  'Slider' },...
  'Style','slider',...
  'Tag','xparamslider');

opt.handles.label.xparam = uicontrol(...
  'Parent',opt.handles.panel.controls,...
  'Units','normalized',...
  'Position',[0.0 0.5 0.625 0.25],...
  'String','xparamLabel',...
  'Style','text',...
  'Tag','xparamLabel');

if ~isempty(opt.ydim)
  opt.handles.slider.yparam = uicontrol(...
    'Parent',opt.handles.panel.controls,...
    'Units','normalized',...
    'BackgroundColor',[0.9 0.9 0.9],...
    'Callback',@cb_slider,...
    'Position',[0.0 0.25 0.625 0.25],...
    'SliderStep', [1/numel(opt.yvalues) 1/numel(opt.yvalues)],...  
    'String',{  'Slider' },...
    'Style','slider',...
    'Tag','yparamSlider');

  opt.handles.label.yparam = uicontrol(...
    'Parent',opt.handles.panel.controls,...
    'Units','normalized',...
    'Position',[0.0 0.0 0.625 0.25],...
    'String','yparamLabel',...
    'Style','text',...
    'Tag','text3');
end

opt.handles.axes.colorbar = axes(...
  'Parent',opt.handles.panel.settings,...
  'Units','normalized',...
  'Position',[0 0 1.0 1],...
  'ButtonDownFcn',@cb_color,...
  'HandleVisibility', 'on', ...
  'Tag','colorbarAxes', ...
  'YLim', [0 1], ...
  'YLimMode', 'manual', ...
  'XLim', [0 1], ...
  'XLimMode', 'manual');

opt.handles.colorbar = colorbar('ButtonDownFcn', []);
set(opt.handles.colorbar, 'Parent', opt.handles.panel.settings);
set(opt.handles.colorbar, 'Position', [0.4 0.2 0.2 0.65]);
nColors = size(colormap, 1);
set(opt.handles.colorbar, 'YTick', [1 nColors/4 nColors/2 3*nColors/4 nColors]);
if (gcf~=opt.handles.figure)
  close gcf; % sometimes there is a new window that opens up
end
%
% set lines
%YLim = get(opt.handles.colorbar, 'YLim');
opt.handles.lines.upperColor = line([-1 0], [32 32], ...
  'Color', 'black', ...
  'LineWidth', 4, ...
  'ButtonDownFcn', @cb_startDrag, ...
  'Parent', opt.handles.colorbar, ...
  'Visible', 'off', ...
  'Tag', 'upperColor');

opt.handles.lines.lowerColor = line([1 2], [32 32], ...
  'Color', 'black', ...
  'LineWidth', 4, ...
  'ButtonDownFcn', @cb_startDrag, ...
  'Parent', opt.handles.colorbar, ...
  'Visible', 'off', ...
  'Tag', 'lowerColor');

opt.handles.menu.colormap = uicontrol(...
  'Parent',opt.handles.panel.settings,...
  'Units','normalized',...
  'BackgroundColor',[1 1 1],...
  'Callback',@cb_colormap,...
  'Position',[0.25 0.95 0.5 0.05],...
  'String',{  'jet'; 'hot'; 'cool'; 'ikelvin'; 'ikelvinr' },...
  'Style','popupmenu',...
  'Value',1, ...
  'Tag','colormapMenu');

opt.handles.checkbox.automatic = uicontrol(...
  'Parent',opt.handles.panel.settings,...
  'Units','normalized',...
  'Callback',@cb_colorbar,...
  'Position',[0.10 0.05 0.75 0.05],...
  'String','automatic',...
  'Style','checkbox',...
  'Value', 1, ... % TODO make this dependet on cfg.zlim
  'Tag','autoCheck');

opt.handles.checkbox.symmetric = uicontrol(...
  'Parent',opt.handles.panel.settings,...
  'Units','normalized',...
  'Callback',@cb_colorbar,...
  'Position',[0.10 0.0 0.75 0.05],...
  'String','symmetric',...
  'Style','checkbox',...
  'Value', 0, ... % TODO make this dependet on cfg.zlim
  'Tag','symCheck');

% buttons

opt.handles.button.decrLowerColor = uicontrol(...
  'Parent', opt.handles.panel.settings,...
  'Units','normalized',...
  'Enable', 'off', ...
  'Callback',@cb_colorbar,...
  'Position',[0.35 0.125 0.15 0.05],...
  'String','-',...
  'Tag','decrLowerColor');

opt.handles.button.incrLowerColor = uicontrol(...
  'Parent', opt.handles.panel.settings,...
  'Units','normalized',...
  'Enable', 'off', ...
  'Callback',@cb_colorbar,...
  'Position',[0.5 0.125 0.15 0.05],...
  'String','+',...
  'Tag','incrLowerColor');


opt.handles.button.decrUpperColor = uicontrol(...
  'Parent', opt.handles.panel.settings,...
  'Units','normalized',...
  'Enable', 'off', ...
  'Callback',@cb_colorbar,...
  'Position',[0.35 0.875 0.15 0.05],...
  'String','-',...
  'Tag','decrUpperColor');

opt.handles.button.incrUpperColor = uicontrol(...
  'Parent', opt.handles.panel.settings,...
  'Units','normalized',...
  'Enable', 'off', ...
  'Callback',@cb_colorbar,...
  'Position',[0.5 0.875 0.15 0.05],...
  'String','+',...
  'Tag','incrUpperColor');

% Handle to the axes that will contain the geometry
opt.handles.axes.movie = axes(...
  'Parent',             opt.handles.panel.visualization_geometry,...
  'Position',           [0 0 1 1],                     ...
  'CameraPosition',     [0.5 0.5 9.16025403784439],    ...
  'CameraPositionMode', get(0,'defaultaxesCameraPositionMode'),...
  'CLim',               get(0,'defaultaxesCLim'),      ...
  'CLimMode',           'manual',                      ...
  'Color',              [0.9 0.9 0.94],                ...
  'ColorOrder',         get(0,'defaultaxesColorOrder'),...
  'XColor',             get(0,'defaultaxesXColor'),    ...
  'YColor',             get(0,'defaultaxesYColor'),    ...
  'ZColor',             get(0,'defaultaxesZColor'),    ...
  'HandleVisibility',   'on',                          ...
  'ButtonDownFcn',      @cb_view,                      ...
  'Tag',                'geometry');

% Handle to the axes that will contain the non-geometry
opt.handles.axes.other = axes(...
  'Parent',             opt.handles.panel.visualization,...
  'Position',           [0 0 1 1],                     ...
  'CameraPosition',     [0.5 0.5 9.16025403784439],    ...
  'CameraPositionMode', get(0,'defaultaxesCameraPositionMode'),...
  'CLim',               get(0,'defaultaxesCLim'),      ...
  'CLimMode',           'manual',                      ...
  'Color',              [0.9 0.9 0.94],                ...
  'ColorOrder',         get(0,'defaultaxesColorOrder'),...
  'XColor',             get(0,'defaultaxesXColor'),    ...
  'YColor',             get(0,'defaultaxesYColor'),    ...
  'ZColor',             get(0,'defaultaxesZColor'),    ...
  'HandleVisibility',   'on',                          ...
  'ButtonDownFcn',      @cb_view,                      ...
  'Tag',                'other');


% Disable axis labels
axis(opt.handles.axes.movie,    'equal');
axis(opt.handles.axes.colorbar, 'equal');
% axis(opt.handles.axes.A, 'equal');
% axis(opt.handles.axes.B, 'equal');
% axis(opt.handles.axes.C, 'equal');
%
axis(opt.handles.axes.movie,    'off');
axis(opt.handles.axes.colorbar, 'off');
% axis(opt.handles.axes.A, 'off');
% axis(opt.handles.axes.B, 'off');
% axis(opt.handles.axes.C, 'off');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_panels(opt)
  for i=1:numel(opt.dat)
    if ~any(opt.issource)
      set(opt.handles.grid{i}, 'cdata', griddata(opt.chanx{i}, opt.chany{i}, opt.mask{i}(:,opt.valx,opt.valy).*opt.dat{i}(:,opt.valx,opt.valy), opt.xdata{i}, opt.nanmask{i}.*opt.ydata{i}, 'v4'));
    else
      set(opt.handles.mesh{i}, 'FaceVertexCData',     squeeze(opt.dat{i}(:,opt.valx,opt.valy)));
      set(opt.handles.mesh{i}, 'FaceVertexAlphaData', squeeze(opt.mask{i}(:,opt.valx,opt.valy)));
    end
  end

  if opt.doplot
    opt.valz
    set(get(opt.handles.axes.other,'children'), 'ydata', opt.dat{1}(opt.valz,:));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACK FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_speed(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);

val = get(h, 'String');
% val = get(opt.handles.slider.speed, 'value');
% val = exp(log(opt.MAX_SPEED) * (val-.5)./0.5);
% 
% % make sure we can easily get back to normal speed
% if abs(val-opt.AVG_SPEED) < 0.08
%   val = opt.AVG_SPEED;
%   set(opt.handles.slider.speed, 'value', 0.5);
% end

speed = str2num(val);
if isempty(speed)
  speed = opt.speed;
end
opt.speed = speed;

set(h, 'String', opt.speed)

% if val >=100
%   set(opt.handles.label.speed, 'String', ['Speed ' num2str(opt.speed, '%.1f'), 'x'])
% elseif val >= 10
%   set(opt.handles.label.speed, 'String', ['Speed ' num2str(opt.speed, '%.2f'), 'x'])
% else
%   set(opt.handles.label.speed, 'String', ['Speed ' num2str(opt.speed, '%.3f'), 'x'])
% end

guidata(h, opt);
uiresume;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACK FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_timer(obj, info, h)
if ~ishandle(h)
  return
end
opt = guidata(h);
delta = opt.speed/numel(opt.xvalues);
val = get(opt.handles.slider.xparam, 'value');
val = val + delta;

if opt.record
  if val>1
    % stop recording
    stop(opt.timer);
    % reset again
    val = 0;    
    set(opt.handles.slider.xparam, 'value', val);
    cb_slider(h);
    cb_recordbutton(opt.handles.button.record);
    % TODO FIXME add some message here
    guidata(h, opt);
    return;
  end
end

while val>1
    val = val-1;  
end
set(opt.handles.slider.xparam, 'value', val);
cb_slider(h);

if opt.record
  pause(.1);
  drawnow;
  vs = version('-release');
  vs = vs(1:4);
  % get starting position via parameter panel    
  currFrame = getframe(opt.handles.figure);
  for i=1:opt.samperframe
    if vs<2010
      opt.vidObj = addframe(opt.vidObj, currFrame);
    else
      writeVideo(opt.vidObj, currFrame);
    end
    
  end
  guidata(h, opt);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACK FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_slider(h, eventdata)
opt = guidata(h);
valx = get(opt.handles.slider.xparam, 'value');
valx = round(valx*(numel(opt.xvalues)-1))+1;
valx = min(valx, numel(opt.xvalues));
valx = max(valx, 1);
set(opt.handles.label.xparam, 'String', [opt.xparam ' ' num2str(opt.xvalues(valx), '%.2f') 's']);


if ~isempty(opt.yvalues)
  valy = get(opt.handles.slider.yparam, 'value');
  valy = round(valy*(numel(opt.yvalues)-1))+1;
  valy = min(valy, numel(opt.yvalues));
  valy = max(valy, 1);
  
  if valy ~= opt.valy
    cb_colorbar(h);
  end
  
  set(opt.handles.label.yparam, 'String', [opt.yparam ' ' num2str(opt.yvalues(valy), '%.2f') 'Hz']);
else
  valy = 1;
end

opt.valx = valx;
opt.valy = valy;

guidata(h, opt);

update_panels(opt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACK FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_playbutton(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
switch get(h, 'string')
  case 'Play'
    set(h, 'string', 'Stop');
    start(opt.timer);
  case 'Stop'
    set(h, 'string', 'Play');
    stop(opt.timer);
end
uiresume;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACK FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_recordbutton(h, eventdata)
if ~ishandle(h)
  return
end

opt = guidata(h);

 
opt.record = ~opt.record; % switch state
guidata(h, opt);

if opt.record
  
  if ~opt.fixedframesfile
    % open a save-file dialog
    [FileName,PathName,FilterIndex] = uiputfile('*.avi', 'Save AVI-file' , 'ft_movie');
    
    if (FileName == 0 & PathName == 0) % aborted
      cb_recordbutton(h);
      return;
    end
    
    opt.framesfile = fullfile(PathName, FileName);
    
    if (FilterIndex==1) % remove .avi again (4 chars)
      opt.framesfile = opt.framesfile(1:end-4);
    end
    
  end
  
  % FIXME open new window to play in there, so that frame getting works
  vs = version('-release');
  vs = vs(1:4);
  if vs<2010
    opt.vidObj = avifile(opt.framesfile, 'FPS', opt.framerate);
  else
    opt.vidObj = VideoWriter(opt.framesfile, 'Uncompressed AVI');
    opt.vidObj.FrameRate = opt.framerate;
    open(opt.vidObj);
  end
 
  %set(opt.handles.figure,'renderer','opengl')
  %opengl software;
  set(opt.handles.figure,'renderer','zbuffer');
  %opt.vidObj = avifile(opt.framesfile, 'fps', opt.framerate, 'quality', 75);
  
  set(h, 'string', 'Stop');  
  guidata(h, opt);
  start(opt.timer);
else
  % FIXME set handle back to old window
  stop(opt.timer);
  if ~isempty(opt.framesfile)
    vs = version('-release');
    vs = vs(1:4);
    if vs<2010
      opt.vidObj = close(opt.vidObj); 
    else    
      close(opt.vidObj);
    end
  end
  set(h, 'string', 'Record');  
  guidata(h, opt);
    
  if (opt.quit)
    close(opt.handles.figure);
  end
end


% 
% % This function should open a new window, plot in there, extract every
% % frame, store the movie in opt and return again
% 
% % this is not needed, no new window is needed
% %scrsz = get(0, 'ScreenSize');
% %f = figure('Position',[1 1 scrsz(3) scrsz(4)]);
% 
% % FIXME disable buttons (apart from RECORD) when recording
% % if record is pressed, stop recording and immediately return
% 
% % adapted from ft_movieplotTFR
% 
% % frequency/time selection
% if ~isempty(opt.yvalues) && any(~isnan(yvalues))
%   if ~isempty(cfg.movietime)
%     indx = cfg.movietime;
%     for iFrame = 1:floor(size(opt.dat, opt.xdim)/cfg.samperframe)
%       indy = ((iFrame-1)*cfg.samperframe+1):iFrame*cfg.samperframe;
%       updateMovie(opt, indx, indy);
%       F(iFrame) = getframe;
%     end
%   elseif ~isempty(cfg.moviefreq)
%     indy = cfg.moviefreq;
%     for iFrame = 1:floor(size(opt.dat, opt.ydim)/cfg.samperframe)
%       indx = ((iFrame-1)*cfg.samperframe+1):iFrame*cfg.samperframe;
%       updateMovie(opt, indx, indy);
%       F(iFrame) = getframe;
%     end
%   else
%     error('Either moviefreq or movietime should contain a bin number')
%   end
% else
%   for iFrame = 1:floor(size(opt.dat, opt.xdim)/cfg.samperframe)
%     indx = ((iFrame-1)*cfg.samperframe+1):iFrame*cfg.samperframe;
%     updateMovie(opt, indx, 1);
%     F(iFrame) = getframe;
%   end
% end
% 
% % save movie
% if ~isempty(cfg.framesfile)
%   save(cfg.framesfile, 'F');
% end
% % play movie
% movie(F, cfg.movierpt, cfg.framespersec);

% uiresume;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACK FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_colormap(h, eventdata)
maps = get(h, 'String');
val = get(h, 'Value');

while ~strcmp(get(h, 'Tag'), 'mainFigure')
  h = get(h, 'Parent');
end

opt =  guidata(h);

cmap = feval(maps{val}, size(colormap, 1));
% if strcmp(maps{val}, 'ikelvin')
%   cmap = ikelvin(size(colormap, 1));
% elseif strcmp(maps{val}, 'kelvin')
%   cmap = kelvin(size(colormap, 1));
% else  
%   cmap = colormap(opt.handles.axes.movie, maps{val});
% end

if get(opt.handles.checkbox.automatic, 'Value')
  colormap(opt.handles.axes.movie_subplot{1}, cmap);
end

adjust_colorbar(opt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACK FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_colorbar(h, eventdata)
if strcmp(get(h, 'Tag'), 'mainFigure') % this is the init call
  incr = false;
  decr = false;
else
  incr = strcmp(get(h, 'String'), '+');
  decr = strcmp(get(h, 'String'), '-');
  lower = strfind(get(h, 'Tag'), 'Lower')>0;
  while ~strcmp(get(h, 'Tag'), 'mainFigure')
    h = get(h, 'Parent');
  end
end

opt =  guidata(h);
if (~incr&&~decr&&exist('eventdata', 'var')&&~isempty(eventdata)) % init call
  caxis(opt.handles.axes.movie, eventdata);
end
zmin = inf;
zmax = -inf;
for i=1:numel(opt.dat)
  [tmpmin tmpmax] = caxis(opt.handles.axes.movie_subplot{i});
  zmin = min(tmpmin, zmin);
  zmax = max(tmpmax, zmax);
end
yLim = get(opt.handles.colorbar, 'YLim');

if incr
  yTick = linspace(zmin, zmax, yLim(end));
  if get(opt.handles.checkbox.symmetric, 'Value')
    zmin = zmin - mean(diff(yTick));
    zmax = zmax + mean(diff(yTick));
  elseif (lower)
    zmin = zmin + mean(diff(yTick));
  else
    zmax = zmax + mean(diff(yTick));
  end
elseif decr
  yTick = linspace(zmin, zmax, yLim(end));
  if get(opt.handles.checkbox.symmetric, 'Value')
    zmin = zmin + mean(diff(yTick));
    zmax = zmax - mean(diff(yTick));
  elseif (lower)
    zmin = zmin - mean(diff(yTick));
  else
    zmax = zmax - mean(diff(yTick));
  end
elseif get(opt.handles.checkbox.automatic, 'Value') % if automatic
  set(opt.handles.lines.upperColor, 'Visible', 'off');
  set(opt.handles.lines.lowerColor, 'Visible', 'off');
  set(opt.handles.lines.upperColor, 'YData', [yLim(end)/2 yLim(end)/2]);
  set(opt.handles.lines.lowerColor, 'YData', [yLim(end)/2 yLim(end)/2]);
  set(opt.handles.button.incrLowerColor, 'Enable', 'off');
  set(opt.handles.button.decrUpperColor, 'Enable', 'off');
  set(opt.handles.button.incrUpperColor, 'Enable', 'off');
  set(opt.handles.button.decrLowerColor, 'Enable', 'off');
    
  if get(opt.handles.checkbox.symmetric, 'Value') % maxabs
    zmax = -inf;
    for i=1:numel(opt.dat)
      tmpmax = max(max(abs(opt.dat{i}(:,:,opt.valy))));
      zmax = max(tmpmax, zmax);
    end
    zmin = -zmax;
  else   % maxmin
    zmax = -inf;
    zmin = inf;
    for i=1:numel(opt.dat)
      tmpmin = min(min(opt.dat{i}(:,:,opt.valy)));
      tmpmax = max(max(opt.dat{i}(:,:,opt.valy)));
      zmax = max(tmpmax, zmax);
      zmin = min(tmpmin, zmin);
    end
    
  end
else
  set(opt.handles.lines.upperColor, 'Visible', 'on');
  set(opt.handles.lines.lowerColor, 'Visible', 'on')
  if get(opt.handles.checkbox.symmetric, 'Value') % maxabs
    set(opt.handles.button.incrLowerColor, 'Enable', 'off');
    set(opt.handles.button.decrUpperColor, 'Enable', 'off');
    set(opt.handles.button.incrUpperColor, 'Enable', 'on');
    set(opt.handles.button.decrLowerColor, 'Enable', 'on');
    for i=1:numel(opt.dat)
      [tmpmin tmpmax] = caxis(opt.handles.axes.movie_subplot{i});
      zmax = max(tmpmax, zmax);
    end
    zmin = -zmax;
  else
    set(opt.handles.button.incrLowerColor, 'Enable', 'on');
    set(opt.handles.button.decrUpperColor, 'Enable', 'on');
    set(opt.handles.button.incrUpperColor, 'Enable', 'on');
    set(opt.handles.button.decrLowerColor, 'Enable', 'on');
    for i=1:numel(opt.dat)
      [tmpmin tmpmax] = caxis(opt.handles.axes.movie_subplot{i});
      zmin = min(tmpmin, zmin);
      zmax = max(tmpmax, zmax);
    end
  end
end % incr, decr, automatic, else

maps = get(opt.handles.menu.colormap, 'String');
cmap = feval(maps{get(opt.handles.menu.colormap, 'Value')}, size(colormap, 1));
for i=1:numel(opt.dat)
  colormap(opt.handles.axes.movie_subplot{i}, cmap);
end

adjust_colorbar(opt);

if (gcf~=opt.handles.figure)
  close gcf; % sometimes there is a new window that opens up
end
for i=1:numel(opt.dat)
  caxis(opt.handles.axes.movie_subplot{i}, [zmin zmax]);
end

nColors = size(colormap, 1);
%yTick = (zmax-zmin)*(get(opt.handles.colorbar, 'YTick')/nColors)+zmin
yTick = (zmax-zmin)*[0 .25 .5 .75 1]+zmin;

%yTick = linspace(zmin, zmax, yLim(end));
% truncate intelligently/-ish
%yTick = get(opt.handles.colorbar, 'YTick')/nColors;
yTick = num2str(yTick', 5);
set(opt.handles.colorbar, 'YTickLabel', yTick, 'FontSize', 8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACK FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_startDrag(h, eventdata)
f = get(h, 'Parent');
while ~strcmp(get(f, 'Tag'), 'mainFigure')
  f = get(f, 'Parent');
end
opt =  guidata(f);
opt.handles.current.line = h;

if strfind(get(h, 'Tag'), 'Color')>0
  opt.handles.current.axes = opt.handles.colorbar;
  opt.handles.current.color = true;
else
  disp('Figure out if it works for xparam and yparam');
  keyboard
end

set(f, 'WindowButtonMotionFcn', @cb_dragLine);

guidata(h, opt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_getposition(h, eventdata)

h   = findobj(h, 'tag', 'mainFigure');
opt = guidata(h);

pos = get(get(h, 'currentaxes'), 'currentpoint');
switch get(get(h, 'currentaxes'), 'tag'),
  case 'geometry'
    if opt.ismesh
      % get the intersection with the mesh
      [ipos, d] = intersect_line(opt.anatomy.pos, opt.anatomy.tri, pos(1,:), pos(2,:));
      [md, ix]  = min(abs(d));
   
      dpos     = opt.anatomy.pos - ipos(ix*ones(size(opt.anatomy.pos,1),1),:);
      opt.valz = nearest(sum(dpos.^2,2),0);
   
    elseif opt.isvolume
    else
    end
    
  case 'other'
  otherwise
end

% if strcmp(get(get(h, 'currentaxes'), 'tag'), 'timecourse')
%   % get the current point
%   pos = get(opt.hy, 'currentpoint');
%   set(opt.sliderx, 'value', nearest(opt.xparam, pos(1,1))./numel(opt.xparam));
%   if isfield(opt, 'hline')
%     set(opt.slidery, 'value', nearest(opt.yparam, pos(1,2))./numel(opt.yparam));
%   end
% elseif strcmp(get(get(h, 'currentaxes'), 'tag'), 'mesh')
%   % get the current point, which is defined as the intersection through the
%   % axis-box (in 3D)
%   pos       = get(opt.hx, 'currentpoint');
%   
%   % get the intersection with the mesh
%   [ipos, d] = intersect_line(opt.pos, opt.tri, pos(1,:), pos(2,:));
%   [md, ix]  = min(abs(d));
%   
%   dpos      = opt.pos - ipos(ix*ones(size(opt.pos,1),1),:);
%   opt.vindx = nearest(sum(dpos.^2,2),0);
%   
%   if isfield(opt, 'parcellation')
%     opt.pindx = find(opt.parcellation(opt.vindx,:));
%     disp(opt.pindx);
%   end
% elseif strcmp(get(get(h, 'currentaxes'), 'tag'), 'mesh2')
%   % get the current point, which is defined as the intersection through the
%   % axis-box (in 3D)
%   pos       = get(opt.hz, 'currentpoint');
%   
%   % get the intersection with the mesh
%   [ipos, d] = intersect_line(opt.pos, opt.tri, pos(1,:), pos(2,:));
%   [md, ix]  = min(abs(d));
%   
%   dpos      = opt.pos - ipos(ix*ones(size(opt.pos,1),1),:);
%   opt.vindx = nearest(sum(dpos.^2,2),0);
%   
%   if isfield(opt, 'parcellation')
%     opt.pindx2 = find(opt.parcellation(opt.vindx,:));
%     disp(opt.pindx2);
%   end
%   
% end
guidata(h, opt);
cb_slider(h);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_dragLine(h, eventdata)
opt =  guidata(h);
pt = get(opt.handles.current.axes, 'CurrentPoint');
yLim = get(opt.handles.colorbar, 'YLim');

% upper (lower) bar must not below (above) lower (upper) bar
if ~(opt.handles.current.line == opt.handles.lines.upperColor && ...
    (any(pt(3)*[1 1]<get(opt.handles.lines.lowerColor, 'YData')) || ...
    yLim(end) <= pt(3))) ...
    && ~(opt.handles.current.line == opt.handles.lines.lowerColor && ...
    (any(pt(3)*[1 1]>get(opt.handles.lines.upperColor, 'YData')) || ...
    yLim(1) >= pt(3)))
  set(opt.handles.current.line, 'YData', pt(3)*[1 1]);
end

adjust_colorbar(opt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_stopDrag(h, eventdata)
while ~strcmp(get(h, 'Tag'), 'mainFigure')
  h = get(h, 'Parent');
end
set(h, 'WindowButtonMotionFcn', '');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function adjust_colorbar(opt)
  % adjust colorbar
  upper = get(opt.handles.lines.upperColor, 'YData');
  lower = get(opt.handles.lines.lowerColor, 'YData');
  if any(round(upper)==0) || any(round(lower)==0)
    return;
  end
  maps = get(opt.handles.menu.colormap, 'String');
  cmap = feval(maps{get(opt.handles.menu.colormap, 'Value')}, size(colormap, 1));
  cmap(round(lower(1)):round(upper(1)), :) = repmat(cmap(round(lower(1)), :), 1+round(upper(1))-round(lower(1)), 1);
  for i=1:numel(opt.dat)
    colormap(opt.handles.axes.movie_subplot{i}, cmap);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = plot_geometry(opt)
  numArgs = numel(opt.dat);
  numRows = floor(sqrt(numArgs));
  numCols = ceil(sqrt(numArgs));
  for i=1:numArgs
    axes(opt.handles.axes.movie);
    opt.handles.axes.movie_subplot{i} = gca;
    if isfield(opt, 'anatomy') && opt.ismesh
      if isfield(opt.anatomy, 'sulc') && ~isempty(opt.anatomy.sulc)
        vdat = opt.anatomy.sulc;
        vdat(vdat>0.5) = 0.5;
        vdat(vdat<-0.5)= -0.5;
        vdat = vdat-min(vdat);
        vdat = 0.35.*(vdat./max(vdat))+0.3;
        vdat = repmat(vdat,[1 3]);
        mesh = ft_plot_mesh(opt.anatomy, 'edgecolor', 'none', 'vertexcolor', vdat);
      else
        mesh = ft_plot_mesh(opt.anatomy, 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5]);
      end
      lighting gouraud
      % set(mesh, 'Parent', opt.handles.axes.movie);
      % mesh = ft_plot_mesh(source, 'edgecolor', 'none', 'vertexcolor', 0*opt.dat(:,1,1), 'facealpha', 0*opt.mask(:,1,1));
      opt.handles.mesh{i} = ft_plot_mesh(opt.anatomy, 'edgecolor', 'none', 'vertexcolor', opt.dat{i}(:,1,1));
      set(opt.handles.mesh{i}, 'AlphaDataMapping', 'scaled');
      set(opt.handles.mesh{i}, 'FaceVertexAlphaData', opt.mask{i}(:,opt.valx,opt.valy));
      % TODO FIXME below does not work
      %set(opt.handles.mesh, 'FaceAlpha', 'flat');
      %set(opt.handles.mesh, 'EdgeAlpha', 'flat');

      lighting gouraud
      cam1 = camlight('left');
%       set(cam1, 'Parent', opt.handles.axes.movie);
      cam2 = camlight('right');
%       set(cam2, 'Parent', opt.handles.axes.movie);
%       set(opt.handles.mesh, 'Parent', opt.handles.axes.movie);
      %   cameratoolbar(opt.handles.figure, 'Show');
    else
      axes(opt.handles.axes.movie)
      [dum, opt.handles.grid{i}] = ft_plot_topo(opt.layout{i}.pos(opt.sellay,1), opt.layout{i}.pos(opt.sellay,2), zeros(numel(opt.sellay{i}),1), 'mask', opt.layout{i}.mask, 'outline', opt.layout{i}.outline, 'interpmethod', 'v4', 'interplim', 'mask', 'parent', opt.handles.axes.movie);
      %[dum, opt.handles.grid] = ft_plot_topo(layout.pos(sellay,1), layout.pos(sellay,2), zeros(numel(sellay),1), 'mask',layout.mask,  'outline', layout.outline, 'interpmethod', 'v4', 'interplim', 'mask', 'parent', opt.handles.axes.movie);
      % set(opt.handles.grid, 'Parent', opt.handles.axes.movie);
      opt.xdata{i}   = get(opt.handles.grid{i}, 'xdata');
      opt.ydata{i}   = get(opt.handles.grid{i}, 'ydata');
      opt.nanmask{i} = 1-get(opt.handles.grid{i}, 'cdata');
      if (gcf~=opt.handles.figure)
        close gcf; % sometimes there is a new window that opens up
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = plot_other(opt)
  dimord  = opt.dimord;
  
  switch dimord
    case {'pos_time' 'pos_freq' 'chan_time' 'chan_freq'}
      opt.doplot    = true;
      opt.doimagesc = false;
    case {'pos_freq_time' 'chan_freq_time' 'chan_chan_freq' 'chan_chan_time' 'pos_pos_freq' 'pos_pos_time'}
      opt.doplot    = false;
      opt.doimagesc = false;
    otherwise
  end
  
  if opt.doplot
    plot(opt.handles.axes.other, opt.xvalues, nanmean(opt.dat{1}(opt.valz,:),1));
  elseif opt.doimagesc
  end
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = ikelvin(m)
%  pos    hue   sat   value
cu = [
  0.0     1/2   0     1.0
  0.125   1/2   0.6   0.95
  0.375   2/3   1.0   0.8
  0.5     2/3   1.0   0.3
  ];

cl = cu;
cl(:, 3:4) = cl(end:-1:1, 3:4);
cl(:, 2)   = cl(:, 2) - 0.5;
cu(:,1)    = cu(:,1)+.5;

x = linspace(0, 1, m)';
l = (x < 0.5); u = ~l;
for i = 1:3
  h(l, i) = interp1(cl(:, 1), cl(:, i+1), x(l));
  h(u, i) = interp1(cu(:, 1), cu(:, i+1), x(u));
end
c = hsv2rgb(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = ikelvinr(m)
%  pos    hue   sat   value
cu = [
  0.0     1/2   0     1.0
  0.125   1/2   0.6   0.95
  0.375   2/3   1.0   0.8
  0.5     2/3   1.0   0.3
  ];

cl = cu;
cl(:, 3:4) = cl(end:-1:1, 3:4);
cl(:, 2)   = cl(:, 2) - 0.5;
cu(:,1)    = cu(:,1)+.5;

x = linspace(0, 1, m)';
l = (x < 0.5); u = ~l;
for i = 1:3
  h(l, i) = interp1(cl(:, 1), cl(:, i+1), x(l));
  h(u, i) = interp1(cu(:, 1), cu(:, i+1), x(u));
end
c = hsv2rgb(h);

c = flipud(c);
end
