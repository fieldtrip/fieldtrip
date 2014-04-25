function moviefunction(cfg, data)
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

ft_defaults
ft_preamble init
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar data

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim',  'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'renamed',    {'zparam', 'funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'parameter', 'funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'mask',      'maskparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'framespersec',      'framerate'});
cfg = ft_checkconfig(cfg, 'deprecated', {'xparam'});
cfg = ft_checkconfig(cfg, 'forbidden',  {'yparam'});

% set data flags
opt = [];
opt.valx = 1;
opt.valy = 1;
opt.issource    = ft_datatype(data, 'source');
opt.isfreq      = ft_datatype(data, 'freq');
opt.istimelock  = ft_datatype(data, 'timelock');

if opt.issource+opt.isfreq+opt.istimelock ~= 1
  error('data cannot be definitely identified as frequency, timelock or source data');
end

if opt.issource
  % transfer anatomical information to opt
  if isfield(data, 'sulc')
    opt.anatomy.sulc = data.sulc;
  end
  
  if isfield(data, 'tri')
    opt.anatomy.tri = data.tri;
  else
    error('source.tri missing, this function requires a triangulated cortical sheet as source model');
  end
  
  if isfield(data, 'pnt')
    opt.anatomy.pnt = data.pnt;
  elseif isfield(data, 'pos')
    opt.anatomy.pos = data.pos;
  else
    error('source.pos or source.pnt is missing');
  end
else
  % identify the interpretation of the functional data
  dtype = ft_datatype(data); 
  switch dtype;
    case 'raw'
      data   = ft_checkdata(data, 'datatype', 'timelock');
      dtype  = ft_datatype(data);
      dimord = data.dimord;
    case  {'timelock' 'freq' 'chan' 'unknown'}
      dimord = data.dimord;
    case 'comp'
      dimord = 'chan_comp';
    otherwise
  end
  dimtok = tokenize(dimord, '_');
end


% set the defaults
cfg.xlim            = ft_getopt(cfg, 'xlim', 'maxmin');
cfg.ylim            = ft_getopt(cfg, 'ylim', 'maxmin');
cfg.zlim            = ft_getopt(cfg, 'zlim', 'maxmin');
cfg.xparam          = ft_getopt(cfg, 'xparam', 'time');
cfg.yparam          = ft_getopt(cfg, 'yparam', []);
if isempty(cfg.yparam) && isfield(data, 'freq')
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
cfg.moviefreq       = ft_getopt(cfg, 'moviefreq', []);
cfg.movietime       = ft_getopt(cfg, 'movietime', []);
cfg.movierpt        = ft_getopt(cfg, 'movierpt', 1);
cfg.interactive     = ft_getopt(cfg, 'interactive', 'yes');
opt.samperframe     = ft_getopt(cfg, 'samperframe',  1);
opt.framerate       = ft_getopt(cfg, 'framerate', 5);
cfg.framesfile      = ft_getopt(cfg, 'framesfile',   []);
opt.record          = ~istrue(cfg.interactive);

% read or create the layout that will be used for plotting:
if ~opt.issource && isfield(cfg, 'layout')
  opt.layout = ft_prepare_layout(cfg);
  
  % Handle the bivariate case

  % Check for bivariate metric with 'chan_chan' in the dimord:
  selchan = strmatch('chan', dimtok);
  isfull  = length(selchan)>1;

  % Check for bivariate metric with a labelcmb field:
  haslabelcmb = isfield(data, 'labelcmb');

  if (isfull || haslabelcmb) && isfield(data, cfg.funparameter)
    % A reference channel is required:
    if ~isfield(cfg, 'refchannel')
      error('no reference channel is specified');
    end

    % check for refchannel being part of selection
    if ~strcmp(cfg.refchannel,'gui')
      if haslabelcmb
        cfg.refchannel = ft_channelselection(cfg.refchannel, unique(data.labelcmb(:)));
      else
        cfg.refchannel = ft_channelselection(cfg.refchannel, data.label);
      end
      if (isfull      && ~any(ismember(data.label, cfg.refchannel))) || ...
          (haslabelcmb && ~any(ismember(data.labelcmb(:), cfg.refchannel)))
        error('cfg.refchannel is a not present in the (selected) channels)')
      end
    end

    if ~isfull,
      % Convert 2-dimensional channel matrix to a single dimension:
      if isempty(cfg.directionality)
        sel1 = find(strcmp(cfg.refchannel, data.labelcmb(:,2)));
        sel2 = find(strcmp(cfg.refchannel, data.labelcmb(:,1)));
      elseif strcmp(cfg.directionality, 'outflow')
        sel1 = [];
        sel2 = find(strcmp(cfg.refchannel, data.labelcmb(:,1)));
      elseif strcmp(cfg.directionality, 'inflow')
        sel1 = find(strcmp(cfg.refchannel, data.labelcmb(:,2)));
        sel2 = [];
      end
      fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.funparameter);
      if length(sel1)+length(sel2)==0
        error('there are no channels selected for plotting: you may need to look at the specification of cfg.directionality');
      end
      data.(cfg.funparameter) = data.(cfg.funparameter)([sel1;sel2],:,:);
      data.label     = [data.labelcmb(sel1,1);data.labelcmb(sel2,2)];
      data           = rmfield(data, 'labelcmb');
    else
      % General case
      sel               = match_str(data.label, cfg.refchannel);
      siz               = [size(data.(cfg.funparameter)) 1];
      if strcmp(cfg.directionality, 'inflow') || isempty(cfg.directionality)
        %the interpretation of 'inflow' and 'outflow' depend on
        %the definition in the bivariate representation of the data
        %in FieldTrip the row index 'causes' the column index channel
        %data.(cfg.funparameter) = reshape(mean(data.(cfg.funparameter)(:,sel,:),2),[siz(1) 1 siz(3:end)]);
        sel1 = 1:siz(1);
        sel2 = sel;
        meandir = 2;
      elseif strcmp(cfg.directionality, 'outflow')
        %data.(cfg.funparameter) = reshape(mean(data.(cfg.funparameter)(sel,:,:),1),[siz(1) 1 siz(3:end)]);
        sel1 = sel;
        sel2 = 1:siz(1);
        meandir = 1;

      elseif strcmp(cfg.directionality, 'inflow-outflow')
        % do the subtraction and recursively call the function again
        tmpcfg = cfg;
        tmpcfg.directionality = 'inflow';
        tmpdata = data;
        tmp     = data.(tmpcfg.funparameter);
        siz     = [size(tmp) 1];
        for k = 1:siz(3)
          for m = 1:siz(4)
            tmp(:,:,k,m) = tmp(:,:,k,m)-tmp(:,:,k,m)';
          end
        end
        tmpdata.(tmpcfg.funparameter) = tmp;
        moviefunction(tmpcfg, tmpdata);
        return;

      elseif strcmp(cfg.directionality, 'outflow-inflow')
        % do the subtraction and recursively call the function again
        tmpcfg = cfg;
        tmpcfg.directionality = 'outflow';
        tmpdata = data;
        tmp     = data.(tmpcfg.funparameter);
        siz     = [size(tmp) 1];
        for k = 1:siz(3)
          for m = 1:siz(4)
            tmp(:,:,k,m) = tmp(:,:,k,m)-tmp(:,:,k,m)';
          end
        end
        tmpdata.(tmpcfg.funparameter) = tmp;
        moviefunction(tmpcfg, tmpdata);
        return;

      end
    end
  end

  % select the channels in the data that match with the layout:
  [opt.seldat, opt.sellay] = match_str(data.label, opt.layout.label);
  
  if isempty(opt.seldat)
    error('labels in data and labels in layout do not match');
  end
  
  
  selcfg = [];
  selcfg.channel = data.label(opt.seldat);
  data = ft_selectdata(selcfg, data);
  
  % get the x and y coordinates and labels of the channels in the data
  opt.chanx = opt.layout.pos(opt.sellay,1);
  opt.chany = opt.layout.pos(opt.sellay,2);
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
opt.framesfile = cfg.framesfile;
opt.fixedframesfile = ~isempty(opt.framesfile);

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
  
  opt.xdim = 2;
  if ~isempty(cfg.yparam)
    opt.ydim = 3;
  else
    opt.ydim = [];
  end

elseif opt.issource
  % TODO FIXME this is hardcoded now, can it be made more generic?
  opt.xdim = 2;
  opt.xvalues = data.time;
  opt.ydim = [];
end

if opt.issource
  if size(data.pos)~=size(opt.dat,1)
    error('inconsistent number of vertices in the cortical mesh');
  end
  
  if ~isfield(data, 'tri')
    error('source.tri missing, this function requires a triangulated cortical sheet as source model');
  end
  
end


if ~isempty(cfg.maskparameter) && ischar(cfg.maskparameter)
  opt.mask = double(getsubfield(data, cfg.maskparameter)~=0);
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

opt = prepareBrainplot(opt);

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
opt.handles.figure = figure(...
  'Units', 'normalized', ...
  'Name','ft_movieplot',...
  'Menu', 'none', ...
  'Toolbar', 'figure', ...
  'NumberTitle','off',...
  'UserData',[],...
  'Tag','mainFigure',...
  'Visible','on');
  %'WindowButtonUpFcn', @cb_stopDrag, ...

%%  panels 
clf


% visualization panel
opt.handles.panel.visualization = uipanel(...
  'tag',    'mainPanels', ... % tag according to position
  'parent', opt.handles.figure,...
  'units',  'normalized', ...
  'title',  'Visualization',...
  'clipping','on',...
  'visible', 'on');

% rearrange 
ft_uilayout(opt.handles.figure, ...
  'tag', 'mainPanels', ...
  'backgroundcolor', [.8 .8 .8], ...
  'hpos', 'auto', ...
  'vpos', .15, ...
  'halign', 'left', ...
  'width', 1, ...
  'height', 0.1);

% settings panel (between different views can be switched)
opt.handles.panel.view = uipanel(...
  'tag',    'sidePanels', ... % tag according to position
  'parent', opt.handles.figure,...
  'units',  'normalized', ...
  'title',  'View',...
  'clipping','on',...
  'visible', 'on');

% settings panel ()
opt.handles.panel.settings = uipanel(...
  'tag',    'sidePanels', ... % tag according to position
  'parent', opt.handles.figure,...
  'units',  'normalized', ...
  'title',  'Settings',...
  'clipping','on',...
  'visible', 'on');

% rearrange panels
ft_uilayout(opt.handles.figure, ...
  'tag', 'sidePanels', ...
  'backgroundcolor', [.8 .8 .8], ...
  'hpos', 'auto', ...
  'vpos', 0.15, ...
  'halign', 'right', ...
  'width', 0.15, ...
  'height', 0.85);

% control panel
opt.handles.panel.controls = uipanel(...
  'tag',    'lowerPanels', ... % tag according to position
  'parent', opt.handles.figure,...
  'units',  'normalized', ...
  'title',  'Controls',...
  'clipping','on',...
  'visible', 'on');

% rearrange 
ft_uilayout(opt.handles.figure, ...
  'tag', 'lowerPanels', ...
  'backgroundcolor', [.8 .8 .8], ...
  'hpos', 'auto', ...
  'vpos', 0, ...
  'halign', 'right', ...
  'width', 1, ...
  'height', 0.15);

ft_uilayout(opt.handles.figure, ...
  'tag', 'mainPanels',  ...
  'retag', 'sidePanels');

ft_uilayout(opt.handles.figure, ...
  'tag', 'sidePanels', ...
  'backgroundcolor', [.8 .8 .8], ...
  'hpos', 'auto', ...
  'vpos', 'align', ...
  'halign', 'left', ...
  'valign', 'top', ...
  'height', .85);

% add axes
% 3 axes for switching to topo viewmode
opt.handles.axes.A = axes(...
  'Parent',opt.handles.panel.view,...
  'Position',[0 0 1 .33], ...
  'Tag','topoAxes',...
  'HandleVisibility', 'on', ...
  'UserData',[]);

opt.handles.axes.B = axes(...
  'Parent',opt.handles.panel.view,...
  'Position',[0 0.33 1 .33], ...
  'Tag','topoAxes',...
  'HandleVisibility', 'on', ...
  'UserData',[]);

opt.handles.axes.C = axes(...
  'Parent',opt.handles.panel.view,...
  'Position',[0 0.66 1 .33], ...
  'Tag','topoAxes',...
  'HandleVisibility', 'on', ...
  'UserData',[]);

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

%
opt.handles.axes.movie = axes(...
  'Parent',opt.handles.panel.visualization,...
  'Position',[0 0 1 1],...
  'CameraPosition',[0.5 0.5 9.16025403784439],...
  'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
  'CLim',get(0,'defaultaxesCLim'),...
  'CLimMode','manual',...
  'Color',[0.9 0.9 0.94],...
  'ColorOrder',get(0,'defaultaxesColorOrder'),...
  'XColor',get(0,'defaultaxesXColor'),...
  'YColor',get(0,'defaultaxesYColor'),...
  'ZColor',get(0,'defaultaxesZColor'),...
  'HandleVisibility', 'on', ...
  'ButtonDownFcn',@cb_view,...
  'Tag','movieAxes');

% Disable axis labels
axis(opt.handles.axes.movie, 'equal');
axis(opt.handles.axes.colorbar, 'equal');
axis(opt.handles.axes.A, 'equal');
axis(opt.handles.axes.B, 'equal');
axis(opt.handles.axes.C, 'equal');
%
axis(opt.handles.axes.movie, 'off');
axis(opt.handles.axes.colorbar, 'off');
axis(opt.handles.axes.A, 'off');
axis(opt.handles.axes.B, 'off');
axis(opt.handles.axes.C, 'off');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateMovie(opt)
if ~opt.issource
  set(opt.handles.grid, 'cdata', griddata(opt.chanx, opt.chany, opt.mask(:,opt.valx,opt.valy).*opt.dat(:,opt.valx,opt.valy), opt.xdata, opt.nanmask.*opt.ydata, 'v4'));
else
  set(opt.handles.mesh, 'FaceVertexCData',  squeeze(opt.dat(:,opt.valx,opt.valy)));
  set(opt.handles.mesh, 'FaceVertexAlphaData', squeeze(opt.mask(:,opt.valx,opt.valy)));
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
delta = opt.speed/size(opt.dat,opt.xdim);
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
valx = round(valx*(size(opt.dat,opt.xdim)-1))+1;
valx = min(valx, size(opt.dat,opt.xdim));
valx = max(valx, 1);
set(opt.handles.label.xparam, 'String', [opt.xparam ' ' num2str(opt.xvalues(valx), '%.2f') 's']);


if ~isempty(opt.yvalues)
  valy = get(opt.handles.slider.yparam, 'value');
  valy = round(valy*(size(opt.dat, opt.ydim)-1))+1;
  valy = min(valy, size(opt.dat, opt.ydim));
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

updateMovie(opt);
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
  colormap(opt.handles.axes.movie, cmap);
else
  adjust_colorbar(opt);
end
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
    zmax = max(max(abs(opt.dat(:,:,opt.valy))));
    zmin = -zmax;
  else   % maxmin
    zmin = min(min(opt.dat(:,:,opt.valy)));
    zmax = max(max(opt.dat(:,:,opt.valy)));
  end
else
  set(opt.handles.lines.upperColor, 'Visible', 'on');
  set(opt.handles.lines.lowerColor, 'Visible', 'on')
  if get(opt.handles.checkbox.symmetric, 'Value') % maxabs
    set(opt.handles.button.incrLowerColor, 'Enable', 'off');
    set(opt.handles.button.decrUpperColor, 'Enable', 'off');
    set(opt.handles.button.incrUpperColor, 'Enable', 'on');
    set(opt.handles.button.decrLowerColor, 'Enable', 'on');
    zmax = max(abs(caxis(opt.handles.axes.movie)));
    zmin = -zmax;
  else
    set(opt.handles.button.incrLowerColor, 'Enable', 'on');
    set(opt.handles.button.decrUpperColor, 'Enable', 'on');
    set(opt.handles.button.incrUpperColor, 'Enable', 'on');
    set(opt.handles.button.decrLowerColor, 'Enable', 'on');
    [zmin zmax] = caxis(opt.handles.axes.movie);
  end
end % incr, decr, automatic, else

maps = get(opt.handles.menu.colormap, 'String');
cmap = feval(maps{get(opt.handles.menu.colormap, 'Value')}, size(colormap, 1));
colormap(opt.handles.axes.movie, cmap);

adjust_colorbar(opt);

if (gcf~=opt.handles.figure)
  close gcf; % sometimes there is a new window that opens up
end

caxis(opt.handles.axes.movie, [zmin zmax]);
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
colormap(opt.handles.axes.movie, cmap);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = prepareBrainplot(opt)
if opt.issource
  if isfield(opt, 'sulc')
    vdat = opt.sulc;
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
  set(mesh, 'Parent', opt.handles.axes.movie);
  % mesh = ft_plot_mesh(source, 'edgecolor', 'none', 'vertexcolor', 0*opt.dat(:,1,1), 'facealpha', 0*opt.mask(:,1,1));
  opt.handles.mesh = ft_plot_mesh(opt.anatomy, 'edgecolor', 'none', 'vertexcolor', opt.dat(:,1,1));
  set(opt.handles.mesh, 'AlphaDataMapping', 'scaled');
  set(opt.handles.mesh, 'FaceVertexAlphaData', opt.mask(:,opt.valx,opt.valy));
  % TODO FIXME below does not work
  %set(opt.handles.mesh, 'FaceAlpha', 'flat');
  %set(opt.handles.mesh, 'EdgeAlpha', 'flat');
  
  lighting gouraud
  cam1 = camlight('left');
  set(cam1, 'Parent', opt.handles.axes.movie);
  cam2 = camlight('right');
  set(cam2, 'Parent', opt.handles.axes.movie);
  set(opt.handles.mesh, 'Parent', opt.handles.axes.movie);
  %   cameratoolbar(opt.handles.figure, 'Show');
else
  axes(opt.handles.axes.movie)
  [dum, opt.handles.grid] = ft_plot_topo(opt.layout.pos(opt.sellay,1), opt.layout.pos(opt.sellay,2), zeros(numel(opt.sellay),1), 'mask', opt.layout.mask, 'outline', opt.layout.outline, 'interpmethod', 'v4', 'interplim', 'mask', 'parent', opt.handles.axes.movie);
  %[dum, opt.handles.grid] = ft_plot_topo(layout.pos(sellay,1), layout.pos(sellay,2), zeros(numel(sellay),1), 'mask',layout.mask,  'outline', layout.outline, 'interpmethod', 'v4', 'interplim', 'mask', 'parent', opt.handles.axes.movie);
  % set(opt.handles.grid, 'Parent', opt.handles.axes.movie);
  opt.xdata   = get(opt.handles.grid, 'xdata');
  opt.ydata   = get(opt.handles.grid, 'ydata');
  opt.nanmask = 1-get(opt.handles.grid, 'cdata');
  if (gcf~=opt.handles.figure)
    close gcf; % sometimes there is a new window that opens up
  end
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
