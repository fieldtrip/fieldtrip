function moviefunction(cfg, data)
% we need cfg.plotfun to plot the data
% data needs to be 3D, N x time x freq (last can be singleton)
%   N needs to correspond to number of vertices (channels, gridpoints, etc)

ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar data

data = ft_checkdata(data, 'datatype', {'timelock', 'freq', 'source'});

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim',  'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'renamed',    {'zparam', 'funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'parameter', 'funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'mask',      'maskparameter'});
cfg = ft_checkconfig(cfg, 'deprecated', {'xparam'});
cfg = ft_checkconfig(cfg, 'forbidden',  {'yparam'});

% set data flags
opt = [];
opt.issource    = ft_datatype(data, 'source');
opt.isfreq      = ft_datatype(data, 'freq');
opt.istimelock  = ft_datatype(data, 'timelock');

if opt.issource+opt.isfreq+opt.istimelock ~= 1
  error('data cannot be definitely deteced as frequency, timelock or source data');
end

% set the defaults
cfg.xlim            = ft_getopt(cfg, 'xlim', 'maxmin');
cfg.ylim            = ft_getopt(cfg, 'ylim', 'maxmin');
cfg.zlim            = ft_getopt(cfg, 'zlim', 'maxmin');
cfg.xparam          = 'time';
if isfield(data, 'freq')
  cfg.yparam        = 'freq';
else
  cfg.yparam        = [];
end
if opt.isfreq
  cfg.funparameter  = ft_getopt(cfg, 'funparameter', 'powspctrm'); % use power as default
elseif opt.istimelock
  cfg.funparameter  = ft_getopt(cfg, 'funparameter', 'avg'); % use power as default
elseif opt.issource
  cfg.funparameter  = ft_getopt(cfg, 'funparameter', 'avg.pow'); % use power as default
end
cfg.maskparameter   = ft_getopt(cfg, 'maskparameter');
cfg.inputfile       = ft_getopt(cfg, 'inputfile',    []);
cfg.samperframe     = ft_getopt(cfg, 'samperframe',  1);
cfg.framespersec    = ft_getopt(cfg, 'framespersec', 5);
cfg.framesfile      = ft_getopt(cfg, 'framesfile',   []);
cfg.moviefreq       = ft_getopt(cfg, 'moviefreq', []);
cfg.movietime       = ft_getopt(cfg, 'movietime', []);
cfg.movierpt        = ft_getopt(cfg, 'movierpt', 1);
cfg.interactive     = ft_getopt(cfg, 'interactive', 'yes');
dointeractive       = istrue(cfg.interactive);

% read or create the layout that will be used for plotting:
if isfield(cfg, 'layout')
  layout = ft_prepare_layout(cfg);
  % select the channels in the data that match with the layout:
  [seldat, sellay] = match_str(data.label, layout.label);
  if isempty(seldat)
    error('labels in data and labels in layout do not match');
  end
  % get the x and y coordinates and labels of the channels in the data
  opt.chanx = layout.pos(sellay,1);
  opt.chany = layout.pos(sellay,2);
else
  if ~opt.issource
    error('you need to specify a layout in case of freq or timelock data');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.xvalues   = data.(cfg.xparam);
opt.yvalues   = [];
if ~isempty(cfg.yparam)
  opt.yvalues = data.(cfg.yparam);
end
opt.dat       = getsubfield(data, cfg.funparameter);

% check consistency of xparam and yparam
% NOTE: i set two different defaults for the 'chan_time' and the 'chan_freq_time' case
if isfield(data,'dimord')
  % get dimord dimensions
  dims = textscan(data.dimord,'%s', 'Delimiter', '_');
  dims = dims{1};
  
  % remove subject dimension
  rpt_dim = strcmp('rpt', dims) | strcmp('subj', dims);
  if sum(rpt_dim)~=0
    dims(rpt_dim) = [];
    opt.dat = squeeze(nanmean(opt.dat, find(rpt_dim)));
  end
  opt.ydim = find(strcmp(cfg.yparam, dims));
  opt.xdim = find(strcmp(cfg.xparam, dims));
  opt.zdim = setdiff(1:ndims(opt.dat), [opt.ydim opt.xdim]);
  if opt.zdim ~=1
    error('input data does not have the correct format of N x time (x freq)');
  end
  % and permute
  opt.dat = permute(opt.dat, [opt.zdim(:)' opt.xdim opt.ydim]);  
     
end

if opt.issource
  if size(data.pos)~=size(opt.dat,1)
    error('inconsistent number of vertices in the cortical mesh');
  end
  
  if ~isfield(data, 'tri')
    error('source.tri missing, this function requires a triangulated cortical sheet as source model');
  end  
  
end

opt.xdim = 2;
if ~isempty(cfg.yparam)
  opt.ydim = 3;
else
  opt.ydim = [];
end

  
if ~isempty(cfg.maskparameter) && ischar(cfg.maskparameter)
  opt.mask = (getsubfield(data, cfg.maskparameter))~=0;
else
  opt.mask =ones(size(opt.dat));
end

if length(opt.xvalues)~=size(opt.dat, opt.xdim)
  error('inconsistent size of "%s" compared to "%s"', cfg.funparameter, cfg.xparam);
end
if ~isempty(opt.ydim) && length(opt.yvalues)~=size(opt.dat,opt.ydim)
  error('inconsistent size of "%s" compared to "%s"', cfg.funparameter, cfg.yparam);
end

% TODO handle colorbar stuff here

% create GUI and plot stuff
opt = createGUI(opt);

if isempty(cfg.yparam)
  set(opt.handles.label.yparam, 'Visible', 'off');
  set(opt.handles.slider.yparam, 'Visible', 'off');
  opt.yparam = '';
else
  opt.yparam = [upper(cfg.yparam(1)) lower(cfg.yparam(2:end))];
  set(opt.handles.label.yparam, 'String', opt.yparam);
end
opt.xparam = [upper(cfg.xparam(1)) lower(cfg.xparam(2:end))];
set(opt.handles.label.xparam, 'String', opt.xparam);

% set timer
opt.timer = timer;
set(opt.timer, 'timerfcn', {@cb_timer, opt.handles.figure}, 'period', 0.1, 'executionmode', 'fixedSpacing');


if opt.issource
  if isfield(data, 'sulc')
    vdat = data.sulc;
    vdat = vdat-min(vdat)+1;
    vdat = vdat./max(vdat);
    vdat = 0.8.*repmat(vdat,[1 3]);
    mesh = ft_plot_mesh(data, 'edgecolor', 'none', 'vertexcolor', vdat);
  else
    mesh = ft_plot_mesh(data, 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5]);
  end
  lighting gouraud
  set(mesh, 'Parent', opt.handles.axes.movie);
  % mesh = ft_plot_mesh(source, 'edgecolor', 'none', 'vertexcolor', 0*opt.dat(:,1,1), 'facealpha', 0*opt.mask(:,1,1));
  opt.handles.mesh = ft_plot_mesh(data, 'edgecolor', 'none', 'vertexcolor', opt.mask(:,1,1).*opt.dat(:,1,1));
  lighting gouraud
  cam1 = camlight('left');
  set(cam1, 'Parent', opt.handles.axes.movie);
  cam2 = camlight('right');
  set(cam2, 'Parent', opt.handles.axes.movie);
  set(opt.handles.mesh, 'Parent', opt.handles.axes.movie);
%   cameratoolbar(opt.handles.figure, 'Show');
else
  axes(opt.handles.axes.movie)
  [dum, opt.handles.grid] = ft_plot_topo(opt.chanx, opt.chany, zeros(numel(opt.chanx),1), 'mask', layout.mask, 'outline', layout.outline, 'interpmethod', 'v4', 'interplim', 'mask', 'parent', opt.handles.axes.movie);
 % set(opt.handles.grid, 'Parent', opt.handles.axes.movie);
 opt.xdata   = get(opt.handles.grid, 'xdata');
 opt.ydata   = get(opt.handles.grid, 'ydata');
 opt.nanmask = 1-get(opt.handles.grid, 'cdata');
  if (gcf~=opt.handles.figure)
    close gcf; % sometimes there is a new window that opens up
  end
end

opt.speed = 1;
guidata(opt.handles.figure, opt);
% init first screen

cb_slider(opt.handles.figure);
%uicontrol(opt.handles.colorbar);
cb_colorbar(opt.handles.figure);

end

%% **************************************************************
%  ********************* CREATE GUI *****************************
%  **************************************************************
function opt = createGUI(opt)
% main figure
opt.handles.figure = figure(...
  'Units','characters',...
  'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
  'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
  'IntegerHandle','off',...
  'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
  'Name','ft_movieplot',...
  'Menu', 'none', ...
  'Toolbar', 'figure', ...
  'NumberTitle','off',...
  'PaperPosition',get(0,'defaultfigurePaperPosition'),...
  'Position',[103.8 14.0769230769231 180.2 50],...
  'HandleVisibility','callback',...
  'WindowButtonUpFcn', @cb_stopDrag, ...
  'UserData',[],...
  'Tag','mainFigure',...
  'Visible','on');

% view panel (between different views can be switched)
opt.handles.panel.view = uipanel(...
  'Parent',opt.handles.figure,...
  'Title','View',...
  'Tag','viewPanel',...
  'Clipping','on',...
  'Visible', 'off', ...
  'Position',[0.842397336293008 -0.00307692307692308 0.155382907880133 1.00153846153846] );

% axes for switching to topo viewmode
opt.handles.axes.topo = axes(...
  'Parent',opt.handles.panel.view,...
  'Position',[0.0441176470588235 0.728706624605678 0.955882352941177 0.205047318611987],...
  'CameraPosition',[0.5 0.5 9.16025403784439],...
  'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
  'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
  'ColorOrder',get(0,'defaultaxesColorOrder'),...
  'LooseInset',[0.115903614457831 0.107232704402516 0.0846987951807229 0.0731132075471698],...
  'XColor',get(0,'defaultaxesXColor'),...
  'YColor',get(0,'defaultaxesYColor'),...
  'ZColor',get(0,'defaultaxesZColor'),...
  'ButtonDownFcn',@(hObject,eventdata)test_movieplot_export('topoAxes_ButtonDownFcn',hObject,eventdata,guidata(hObject)),...
  'Tag','topoAxes',...
  'HandleVisibility', 'on', ...
  'UserData',[]);

% axes for switching to multiplot viewmode
opt.handles.axes.multi = axes(...
  'Parent',opt.handles.panel.view,...
  'Position',[0.0441176470588235 0.413249211356467 0.955882352941177 0.205047318611987],...
  'CameraPosition',[0.5 0.5 9.16025403784439],...
  'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
  'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
  'ColorOrder',get(0,'defaultaxesColorOrder'),...
  'LooseInset',[0.0302261898791314 0.128446186295888 0.0220883695270575 0.0875769452017415],...
  'XColor',get(0,'defaultaxesXColor'),...
  'YColor',get(0,'defaultaxesYColor'),...
  'ZColor',get(0,'defaultaxesZColor'),...
  'ButtonDownFcn',@(hObject,eventdata)test_movieplot_export('multiAxes_ButtonDownFcn',hObject,eventdata,guidata(hObject)),...
  'Tag','multiAxes',...
   'HandleVisibility', 'on', ...
  'UserData',[]);

% axes for switching to singleplot viewmode
opt.handles.axes.single = axes(...
  'Parent',opt.handles.panel.view,...
  'Position',[0.0441176470588235 0.0977917981072555 0.955882352941177 0.205047318611987],...
  'CameraPosition',[0.5 0.5 9.16025403784439],...
  'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
  'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
  'ColorOrder',get(0,'defaultaxesColorOrder'),...
  'LooseInset',[0.0302261898791314 0.128446186295888 0.0220883695270575 0.0875769452017415],...
  'XColor',get(0,'defaultaxesXColor'),...
  'YColor',get(0,'defaultaxesYColor'),...
  'ZColor',get(0,'defaultaxesZColor'),...
  'Tag','singleAxes',...
   'HandleVisibility', 'on', ...
  'UserData',[]);

% main control panel
opt.handles.panel.control = uipanel(...
  'Parent',opt.handles.figure,...
  'Title','Controls',...
  'Tag','controlPanel',...
  'Clipping','on',...
  'Position',[0.660377358490566 -0.00307692307692308 0.177580466148724 0.153846153846154]);

% buttons
opt.handles.button.play = uicontrol(...
  'Parent',opt.handles.panel.control,...
  'Units','normalized',...
  'Callback',@cb_playbutton,...
  'Position',[0.0973193473193473 0.457831325301205 0.314685314685315 0.530120481927711],...
  'String','Play',...
  'Tag','playButton');

opt.handles.button.record = uicontrol(...
  'Parent',opt.handles.panel.control,...
  'Units','normalized',...
  'Callback',@(hObject,eventdata)test_movieplot_export('recordButton_Callback',hObject,eventdata,guidata(hObject)),...
  'Position',[0.551864801864802 0.457831325301205 0.321678321678322 0.530120481927711],...
  'String','Record',...
  'Tag','recordButton');

% speed control
opt.handles.slider.speed = uicontrol(...
  'Parent',opt.handles.panel.control,...
  'Units','normalized',...
  'BackgroundColor',[0.9 0.9 0.9],...
  'Callback',@cb_speed, ...
  'Position',[0.0256410256410256 0.0843373493975904 0.493589743589744 0.204819277108434],...
  'String',{  'Slider' },...
  'Style','slider',...
  'Value',0.5,...
  'Tag','speedSlider');
opt.MIN_SPEED = 0.01;
opt.AVG_SPEED = 1;
opt.MAX_SPEED = 100;

opt.handles.label.speed = uicontrol(...
  'Parent',opt.handles.panel.control,...
  'Units','normalized',...
  'HorizontalAlignment','left',...
  'Position',[0.532051282051282 0.108433734939759 0.448717948717949 0.168674698795181],...
  'String','Speed 1.0x',...
  'Style','text',...
  'Value',0.5,...
  'Tag','speedLabel');

% parameter control
opt.handles.panel.parameter = uipanel(...
  'Parent',opt.handles.figure,...
  'Title','Parameter',...
  'Tag','parameterPanel',...
  'Clipping','on',...
  'Position',[-0.00110987791342952 -0.00307692307692308 0.658157602663707 0.155384615384615]);

opt.handles.slider.xparam = uicontrol(...
  'Parent',opt.handles.panel.parameter,...
  'Units','normalized',...
  'BackgroundColor',[0.9 0.9 0.9],...
  'Callback',@cb_slider,...
  'Position',[0.0186757215619694 0.250000000000001 0.969439728353141 0.214285714285714],...
  'String',{  'Slider' },...
  'Style','slider',...
  'Tag','xparamslider');

opt.handles.label.xparam = uicontrol(...
  'Parent',opt.handles.panel.parameter,...
  'Units','normalized',...
  'Position',[0.0186757215619694 0.0476190476190483 0.967741935483871 0.178571428571429],...
  'String','xparamLabel',...
  'Style','text',...
  'Tag','xparamLabel');

opt.handles.slider.yparam = uicontrol(...
  'Parent',opt.handles.panel.parameter,...
  'Units','normalized',...
  'BackgroundColor',[0.9 0.9 0.9],...
  'Callback',@cb_slider,...
  'Position',[0.0186757215619694 0.690476190476191 0.969439728353141 0.202380952380952],...
  'String',{  'Slider' },...
  'Style','slider',...
  'Tag','yparamSlider');

opt.handles.label.yparam = uicontrol(...
  'Parent',opt.handles.panel.parameter,...
  'Units','normalized',...
  'Position',[0.0186757215619694 0.500000000000001 0.967741935483871 0.178571428571429],...
  'String','yparamLabel',...
  'Style','text',...
  'Tag','text3');

opt.handles.buttongroup.color = uibuttongroup(...
  'Parent',opt.handles.figure,...
  'Title','Colormap',...
  'Tag','colorGroup',...
  'Clipping','on',...
  'Position',[0.725860155382908 0.150769230769231 0.112097669256382 0.847692307692308],...
  'SelectedObject',[],...
  'SelectionChangeFcn',[],...
  'OldSelectedObject',[]);

opt.handles.axes.colorbar = axes(...
  'Parent',opt.handles.buttongroup.color,...
  'Position',[-0.0103092783505155 0.157303370786517 1.04123711340206 0.771535580524345],...
  'CameraPosition',[0.5 0.5 9.16025403784439],...
  'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
  'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
  'ColorOrder',get(0,'defaultaxesColorOrder'),...
  'LooseInset',[6.37540259335975e-007 1.19260520025353e-007 4.65894804899366e-007 8.13139909263772e-008],...
  'XColor',get(0,'defaultaxesXColor'),...
  'YColor',get(0,'defaultaxesYColor'),...
  'ZColor',get(0,'defaultaxesZColor'),...
  'ButtonDownFcn',@(hObject,eventdata)test_movieplot_export('colorbarAxes_ButtonDownFcn',hObject,eventdata,guidata(hObject)),...
  'HandleVisibility', 'on', ...
  'Tag','colorbarAxes', ...
  'YLim', [0 1], ...
  'YLimMode', 'manual', ...
  'XLim', [0 1], ...
  'XLimMode', 'manual');

% set colorbar
opt.handles.colorbar = colorbar('ButtonDownFcn', []);
set(opt.handles.colorbar, 'Parent', opt.handles.figure);
set(opt.handles.colorbar, 'Position', [0.7725 0.305 0.02 0.6]);
if (gcf~=opt.handles.figure)
  close gcf; % sometimes there is a new window that opens up
end

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
  'Parent',opt.handles.buttongroup.color,...
  'Units','normalized',...
  'BackgroundColor',[1 1 1],...
  'Callback',@cb_colormap,...
  'Position',[0.0721649484536082 0.951310861423221 0.865979381443299 0.0374531835205993],...
  'String',{  'jet'; 'hot'; 'cool' },...
  'Style','popupmenu',...
  'Value',1, ...
  'Tag','colormapMenu');

opt.handles.checkbox.auto = uicontrol(...
  'Parent',opt.handles.buttongroup.color,...
  'Units','normalized',...
  'Callback',@cb_colorbar,...
  'Position',[0.103092783505155 0.0898876404494382 0.77319587628866 0.0430711610486891],...
  'String','automatic',...
  'Style','checkbox',...
  'Value', 1, ... % TODO make this dependet on cfg.zlim
  'Tag','autoCheck');

opt.handles.checkbox.symmetric = uicontrol(...
  'Parent',opt.handles.buttongroup.color,...
  'Units','normalized',...
  'Callback',@cb_colorbar,...
  'Position',[0.103092783505155 0.0449438202247191 0.77319587628866 0.0430711610486891],...
  'String','symmetric',...
  'Style','checkbox',...
  'Value', 0, ... % TODO make this dependet on cfg.zlim
  'Tag','symCheck');

% buttons

opt.handles.button.decrUpperColor = uicontrol(...
  'Parent', opt.handles.buttongroup.color,...
  'Units','normalized',...
  'Enable', 'off', ...
  'Callback',@cb_colorbar,...
  'Position',[0.503092783505155 0.9149438202247191 0.17319587628866 0.0330711610486891],...
  'String','-',...
  'Tag','decrUpperColor');

opt.handles.button.incrLowerColor = uicontrol(...
  'Parent', opt.handles.buttongroup.color,...
  'Units','normalized',...
  'Enable', 'off', ...
  'Callback',@cb_colorbar,...
  'Position',[0.303092783505155 0.1449438202247191 0.17319587628866 0.0330711610486891],...
  'String','+',...
  'Tag','incrLowerColor');

opt.handles.button.incrUpperColor = uicontrol(...
  'Parent', opt.handles.buttongroup.color,...
  'Units','normalized',...
  'Enable', 'off', ...
  'Callback',@cb_colorbar,...
  'Position',[0.303092783505155 0.9149438202247191 0.17319587628866 0.0330711610486891],...
  'String','+',...
  'Tag','incrUpperColor');

opt.handles.button.decrLowerColor = uicontrol(...
  'Parent', opt.handles.buttongroup.color,...
  'Units','normalized',...
  'Enable', 'off', ...
  'Callback',@cb_colorbar,...
  'Position',[0.503092783505155 0.1449438202247191 0.17319587628866 0.0330711610486891],...
  'String','-',...
  'Tag','decrLowerColor');

opt.handles.axes.movie = axes(...
  'Parent',opt.handles.figure,...
  'Position',[-0.00110987791342952 0.150769230769231 0.722530521642619 0.847692307692308],...
  'CameraPosition',[0.5 0.5 9.16025403784439],...
  'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
  'CLim',get(0,'defaultaxesCLim'),...
  'CLimMode','manual',...
  'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
  'ColorOrder',get(0,'defaultaxesColorOrder'),...
  'LooseInset',[0.120510948905109 0.11 0.088065693430657 0.075],...
  'XColor',get(0,'defaultaxesXColor'),...
  'YColor',get(0,'defaultaxesYColor'),...
  'ZColor',get(0,'defaultaxesZColor'),...
  'HandleVisibility', 'on', ...
  'ButtonDownFcn',@(hObject,eventdata)test_movieplot_export('movieAxes_ButtonDownFcn',hObject,eventdata,guidata(hObject)),...
  'Tag','movieAxes');

% Disable axis labels
axis(opt.handles.axes.movie, 'equal');
axis(opt.handles.axes.colorbar, 'equal');
axis(opt.handles.axes.topo, 'equal');
axis(opt.handles.axes.multi, 'equal');
axis(opt.handles.axes.single, 'equal');

axis(opt.handles.axes.movie, 'off');
axis(opt.handles.axes.colorbar, 'off');
axis(opt.handles.axes.topo, 'off');
axis(opt.handles.axes.multi, 'off');
axis(opt.handles.axes.single, 'off');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateMovie(opt, valx, valy)
if ~opt.issource
  set(opt.handles.grid, 'cdata',  griddata(opt.chanx, opt.chany, opt.mask(:,valx,valy).*opt.dat(:,valx,valy), opt.xdata, opt.nanmask.*opt.ydata, 'v4'));
else
  set(opt.handles.mesh, 'FaceVertexCData',  squeeze(opt.mask(:,valx,valy).*opt.dat(:,valx,valy)));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_speed(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);

val = get(opt.handles.slider.speed, 'value');
val = exp(log(opt.MAX_SPEED) * (val-.5)./0.5);

% make sure we can easily get back to normal speed
if abs(val-opt.AVG_SPEED) < 0.08
  val = opt.AVG_SPEED;
  set(opt.handles.slider.speed, 'value', 0.5);
end
  
opt.speed = val;

if val >=100
  set(opt.handles.label.speed, 'String', ['Speed ' num2str(opt.speed, '%.1f'), 'x'])
elseif val >= 10
  set(opt.handles.label.speed, 'String', ['Speed ' num2str(opt.speed, '%.2f'), 'x'])  
else
  set(opt.handles.label.speed, 'String', ['Speed ' num2str(opt.speed, '%.3f'), 'x'])
end

guidata(h, opt);
uiresume;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_timer(obj, info, h)
if ~ishandle(h)
  return
end
opt = guidata(h);
delta = opt.speed/size(opt.dat,opt.xdim);
val = get(opt.handles.slider.xparam, 'value');
val = val + delta;
if val>1
  val = val-1;
end
set(opt.handles.slider.xparam, 'value', val);
cb_slider(h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_slider(h, eventdata)
opt = guidata(h);
valx = get(opt.handles.slider.xparam, 'value');
valx = round(valx*(size(opt.dat,opt.xdim)-1))+1;
valx = min(valx, size(opt.dat,opt.xdim));
valx = max(valx, 1);
set(opt.handles.label.xparam, 'String', [opt.xparam ' ' num2str(opt.xvalues(valx), '%.2f') 's']);


if ~isempty(opt.yvalues)
  valy = get(opt.handles.slider.yparam, 'value');
  valy = round(valy*(size(opt.dat, opt.ydim)-1))+1;
  valy = min(valy, size(opt.dat, opt.ydim));
  valy = max(valy, 1);
  set(opt.handles.label.yparam, 'String', [opt.yparam ' ' num2str(opt.yvalues(valy), '%.2f') 'Hz']);
else
  valy = 1;
end

updateMovie(opt, valx, valy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
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
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_colormap(h, eventdata)
maps = get(h, 'String');
val = get(h, 'Value');

while ~strcmp(get(h, 'Tag'), 'mainFigure')
  h = get(h, 'Parent');
end
  
opt =  guidata(h);

cmap = colormap(opt.handles.axes.movie, maps{val});

if get(opt.handles.checkbox.auto, 'Value')
  colormap(opt.handles.axes.movie, cmap);
else
  adjust_colorbar(opt);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
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
[zmin zmax] = caxis(opt.handles.axes.movie);
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
elseif get(opt.handles.checkbox.auto, 'Value') % if automatic
  set(opt.handles.lines.upperColor, 'Visible', 'off');
  set(opt.handles.lines.lowerColor, 'Visible', 'off');
  set(opt.handles.lines.upperColor, 'YData', [yLim(end)/2 yLim(end)/2]);
  set(opt.handles.lines.lowerColor, 'YData', [yLim(end)/2 yLim(end)/2]);
  set(opt.handles.button.incrUpperColor, 'Enable', 'off');
  set(opt.handles.button.decrUpperColor, 'Enable', 'off');
  set(opt.handles.button.incrLowerColor, 'Enable', 'off');
  set(opt.handles.button.decrLowerColor, 'Enable', 'off');
  if get(opt.handles.checkbox.symmetric, 'Value') % maxabs
    zmax = max(abs(opt.dat(:)));
    zmin = -zmax;
  else % maxmin
    zmin = min(opt.dat(:));
    zmax = max(opt.dat(:));
  end
else
  set(opt.handles.lines.upperColor, 'Visible', 'on');
  set(opt.handles.lines.lowerColor, 'Visible', 'on')
  set(opt.handles.button.incrUpperColor, 'Enable', 'on');
  set(opt.handles.button.decrLowerColor, 'Enable', 'on');
  if get(opt.handles.checkbox.symmetric, 'Value') % maxabs
    set(opt.handles.button.decrUpperColor, 'Enable', 'off');
    set(opt.handles.button.incrLowerColor, 'Enable', 'off');
    zmax = max(abs(caxis(opt.handles.axes.movie)));
    zmin = -zmax;
  else
    set(opt.handles.button.decrUpperColor, 'Enable', 'on');
    set(opt.handles.button.incrLowerColor, 'Enable', 'on');
    [zmin zmax] = caxis(opt.handles.axes.movie);
  end
end % incr, decr, automatic, else

maps = get(opt.handles.menu.colormap, 'String');
cmap = colormap(opt.handles.axes.movie, maps{get(opt.handles.menu.colormap, 'Value')});

adjust_colorbar(opt);

if (gcf~=opt.handles.figure)
  close gcf; % sometimes there is a new window that opens up
end

caxis(opt.handles.axes.movie, [zmin zmax]);  
yTick = linspace(zmin, zmax, yLim(end));
% truncate intelligently/-ish
yTick = yTick(get(opt.handles.colorbar, 'YTick'));
yTick = num2str(yTick', 5);
set(opt.handles.colorbar, 'YTickLabel', yTick, 'FontSize', 8);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
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


function adjust_colorbar(opt)
% adjust colorbar
upper = get(opt.handles.lines.upperColor, 'YData');
lower = get(opt.handles.lines.lowerColor, 'YData');
maps = get(opt.handles.menu.colormap, 'String');
cmap = colormap(opt.handles.axes.movie, maps{get(opt.handles.menu.colormap, 'Value')});
cmap(round(lower(1)):round(upper(1)), :) = repmat(cmap(round(lower(1)), :), 1+round(upper(1))-round(lower(1)), 1);
colormap(opt.handles.axes.movie, cmap);
end

function cb_stopDrag(h, eventdata)
while ~strcmp(get(h, 'Tag'), 'mainFigure')
  h = get(h, 'Parent');
end
set(h, 'WindowButtonMotionFcn', '');
end