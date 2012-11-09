function [cfg, M] = ft_sourcemovie(cfg, source)

% FT_SOURCEMOVIE displays the source reconstruction on a cortical mesh
% and allows the user to scroll through time with a movie
%
% Use as
%   ft_sourcemovie(cfg, source)
% where the input source data is obtained from FT_SOURCEANALYSIS and cfg is
% a configuratioun structure that should contain
%
%  cfg.funparameter    = string, functional parameter that is color coded (default = 'avg.pow')
%  cfg.maskparameter   = string, functional parameter that is used for opacity (default = [])
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SOURCEPLOT, FT_SOURCEINTERPOLATE

% Copyright (C) 2011, Robert Oostenveld
%
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig
ft_preamble loadvar source

% ensure that the input data is valiud for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
source = ft_checkdata(source, 'datatype', 'source', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',	 {'zparam',    'funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	 {'parameter', 'funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	 {'mask',      'maskparameter'});

% these are not needed any more, once the source structure has a proper dimord
% cfg = ft_checkconfig(cfg, 'deprecated', 'xparam');
% cfg = ft_checkconfig(cfg, 'deprecated', 'yparam');

% get the options
xlim          = ft_getopt(cfg, 'xlim');
ylim          = ft_getopt(cfg, 'ylim');
zlim          = ft_getopt(cfg, 'zlim');
olim          = ft_getopt(cfg, 'alim');                           % don't use alim as variable name
xparam        = ft_getopt(cfg, 'xparam', 'time');                 % use time as default
yparam        = ft_getopt(cfg, 'yparam');                         % default is dealt with below
funparameter  = ft_getopt(cfg, 'funparameter', 'avg.pow');        % use power as default
maskparameter = ft_getopt(cfg, 'maskparameter');

if isempty(yparam) && isfield(source, 'freq')
  % the default is freq (if present)
  yparam = 'freq';
end

% update the configuration
cfg.funparameter  = funparameter;
cfg.maskparameter = maskparameter;
cfg.xparam        = xparam;
cfg.yparam        = yparam;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fun = getsubfield(source, funparameter); % might be avg.pow
if size(source.pos)~=size(fun,1)
  error('inconsistent number of vertices in the cortical mesh');
end

if ~isfield(source, 'tri')
  error('source.tri missing, this function requires a triangulated cortical sheet as source model');
end

if ~isempty(maskparameter) && ischar(maskparameter)
  mask = double(getsubfield(source, maskparameter));
else
  mask = 0.5*ones(size(fun));
end

xparam = source.(xparam);
if length(xparam)~=size(fun,2)
  error('inconsistent size of "%s" compared to "%s"', cfg.funparameter, cfg.xparam);
end

if ~isempty(yparam)
  yparam = source.(yparam);
  if length(yparam)~=size(fun,3)
    error('inconsistent size of "%s" compared to "%s"', cfg.funparameter, cfg.yparam);
  end
end

if isempty(xlim)
  xlim(1) = min(xparam);
  xlim(2) = max(xparam);
end

xbeg = nearest(xparam, xlim(1));
xend = nearest(xparam, xlim(2));
% update the configuration
cfg.xlim = xparam([xbeg xend]);

if ~isempty(yparam)
  if isempty(ylim)
    ylim(1) = min(yparam);
    ylim(2) = max(yparam);
  end
  ybeg = nearest(yparam, ylim(1));
  yend = nearest(yparam, ylim(2));
  % update the configuration
  cfg.ylim = xparam([ybeg yend]);
  hasyparam = true;
else
  % this allows us not to worry about the yparam any more
  yparam = nan;
  ybeg = 1;
  yend = 1;
  cfg.ylim = [];
  hasyparam = false;
end

% make a subselection of the data
xparam  = xparam(xbeg:xend);
yparam  = yparam(ybeg:yend);
fun     = fun(:,xbeg:xend,ybeg:yend);
mask    = mask(:,xbeg:xend,ybeg:yend);
clear xbeg xend ybeg yend

if isempty(zlim)
  zlim(1) = min(fun(:));
  zlim(2) = max(fun(:));
  % update the configuration
  cfg.zlim = zlim;
end

if isempty(olim)
  olim(1) = min(mask(:));
  olim(2) = max(mask(:));
  if olim(1)==olim(2)
    olim(1) = 0;
    olim(2) = 1;
  end
  % update the configuration
  cfg.alim = olim;
end

% collect the data and the options to be used in the figure
opt.cfg     = cfg;
opt.xparam  = xparam;
opt.yparam  = yparam;
opt.xval    = 0;
opt.yval    = 0;
opt.dat     = fun;
opt.mask    = mask;
opt.pos     = source.pos;
opt.vindx   = source.inside(:);
opt.speed   = 1;
opt.record  = 0;
opt.threshold = 0;
opt.frame   = 0;
opt.cleanup = false;

% get a handle to a figure
h  = gcf;
set(h, 'color', [1 1 1]);
set(h, 'toolbar', 'figure');
set(h, 'visible', 'on');
set(h, 'CloseRequestFcn', @cb_quitbutton);
set(h, 'position', [100 200 700 500]);
set(h, 'windowbuttondownfcn', @cb_getposition); 

% get timer object
t = timer;
set(t, 'timerfcn', {@cb_timer, h}, 'period', 0.1, 'executionmode', 'fixedSpacing');

% make the user interface elements
cambutton    = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'light', 'userdata', 'C');
playbutton   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'play',   'userdata', 'p');
recordbutton = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'record', 'userdata', 'r');
quitbutton   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'quit',   'userdata', 'q');

thrmin   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'downarrow');
thr      = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'threshold', 'userdata', 't');
thrplus  = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'uparrow');
spdmin   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'shift+downarrow');
spd      = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'speed','userdata', 's');
spdplus  = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'shift+uparrow');
climmin  = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'leftarrow');
clim     = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'colorlim', 'userdata', 'z');
climplus = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'rightarrow');
sliderx  = uicontrol('parent', h, 'units', 'normalized', 'style', 'slider',     'string', sprintf('%s = ', cfg.xparam));
stringx  = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');
slidery  = uicontrol('parent', h, 'units', 'normalized', 'style', 'slider',     'string', sprintf('%s = ', cfg.yparam));
stringy  = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');
stringz  = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');

set(cambutton,    'position', [0.095 0.28 0.09 0.05], 'callback', @cb_keyboard);
set(quitbutton,   'position', [0.005 0.28 0.09 0.05], 'callback', @cb_keyboard);
set(playbutton,   'position', [0.005 0.22 0.09 0.05], 'callback', @cb_keyboard);
set(recordbutton, 'position', [0.095 0.22 0.09 0.05], 'callback', @cb_keyboard);
set(thrmin,       'position', [0.005 0.16 0.03 0.05], 'callback', @cb_keyboard);
set(thr,          'position', [0.035 0.16 0.12 0.05], 'callback', @cb_keyboard);
set(thrplus,      'position', [0.155 0.16 0.03 0.05], 'callback', @cb_keyboard);
set(climmin,      'position', [0.005 0.10 0.03 0.05], 'callback', @cb_keyboard);
set(clim,         'position', [0.035 0.10 0.12 0.05], 'callback', @cb_keyboard);
set(climplus,     'position', [0.155 0.10 0.03 0.05], 'callback', @cb_keyboard);
set(spdmin,       'position', [0.005 0.04 0.03 0.05], 'callback', @cb_keyboard);
set(spd,          'position', [0.035 0.04 0.12 0.05], 'callback', @cb_keyboard);
set(spdplus,      'position', [0.155 0.04 0.03 0.05], 'callback', @cb_keyboard);
set(sliderx,      'position', [0.01 0.4 0.26 0.03], 'callback',  @cb_slider);%[0.200 0.04  0.78 0.03], 'callback', @cb_slider);
set(slidery,      'position', [0.200 0.005 0.78 0.03], 'callback', @cb_slider);
set(stringx,      'position', [0.800 0.93 0.18 0.03]);
set(stringy,      'position', [0.800 0.90 0.18 0.03]);
set(stringz,      'position', [0.650 0.96 0.33 0.03]);

set(stringx, 'string', sprintf('%s = ', cfg.xparam));
set(stringy, 'string', sprintf('%s = ', cfg.yparam));
set(stringz, 'string', sprintf('position = '));
set(stringx, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);
set(stringy, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);
set(stringz, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);

% create axes object to contain the mesh
hx = axes;
set(hx, 'position', [0.4 0.08 0.6 0.8]);
set(hx, 'tag', 'mesh');
if isfield(source, 'sulc')
  vdat = source.sulc;
  vdat = vdat-min(vdat)+1;
  vdat = vdat./max(vdat);
  vdat = 0.8.*repmat(vdat,[1 3]);
  hs1 = ft_plot_mesh(source, 'edgecolor', 'none', 'vertexcolor', vdat);
else
  hs1 = ft_plot_mesh(source, 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5]);
end
lighting gouraud
hs = ft_plot_mesh(source, 'edgecolor', 'none', 'vertexcolor', 0*opt.dat(:,1,1), 'facealpha', 0*opt.mask(:,1,1));
lighting gouraud
cam1 = camlight('left');
cam2 = camlight('right');
caxis(cfg.zlim);
alim(cfg.alim);

% create axis object to contain a time course
hy = axes;
set(hy, 'position', [0.02 0.5 0.3 0.35]);
set(hy, 'yaxislocation', 'right');

if ~hasyparam
  tline = plot(opt.xparam, mean(opt.dat(opt.vindx,:))); hold on;
  abc = axis;
  axis([opt.xparam(1) opt.xparam(end) abc(3:4)]);
  vline = plot(opt.xparam(1)*[1 1], abc(3:4), 'r');
else
  error('not yet implemented');  
end
set(hy, 'tag', 'timecourse');

% remember the various handles
opt.h   = h;  % handle to the figure
opt.hs  = hs; % handle to the mesh
opt.hx  = hx; % handle to the axes containing the mesh
opt.hy  = hy; % handle to the axes containing the timecourse
opt.cam = [cam1 cam2]; % handles to the light objects
opt.vline = vline; % handle to the line in the ERF plot
opt.tline = tline; % handle to the ERF
opt.playbutton   = playbutton; % handle to the playbutton
opt.recordbutton = recordbutton; % handle to the recordbutton
opt.quitbutton   = quitbutton; % handle to the quitbutton

%opt.p   = p;
opt.t   = t;
%opt.hx  = hx;
%opt.hy  = hy;
opt.sliderx  = sliderx;
opt.slidery  = slidery;
opt.stringx  = stringx;
opt.stringy  = stringy;
opt.stringz  = stringz;

if ~hasyparam
  set(opt.slidery, 'visible', 'off');
  set(opt.stringy, 'visible', 'off');
end

setappdata(h, 'opt', opt);

while opt.cleanup==0
  uiwait(h);
  opt = getappdata(h, 'opt');
end
stop(opt.t);

if nargout
  M = opt.movie;
end

delete(h);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_slider(h, eventdata)

persistent previous_valx previous_valy previous_vindx

if isempty(previous_valx)
  previous_valx = 0;
end
if isempty(previous_valy)
  previous_valy = 0;
end

h    = getparent(h);
opt  = getappdata(h, 'opt');
valx = get(opt.sliderx, 'value');
valx = round(valx*(size(opt.dat,2)-1))+1;
valx = min(valx, size(opt.dat,2));
valx = max(valx, 1);

valy = get(opt.slidery, 'value');
valy = round(valy*(size(opt.dat,3)-1))+1;
valy = min(valy, size(opt.dat,3));
valy = max(valy, 1);

mask = squeeze(opt.mask(:,valx,valy));
mask(opt.dat(:,valx,valy)<opt.threshold) = 0;

% update stuff
if previous_valx~=valx || previous_valy~=valy
  % update strings
  set(opt.stringx, 'string', sprintf('%s = %3.3f\n', opt.cfg.xparam, opt.xparam(valx)));
  set(opt.stringy, 'string', sprintf('%s = %3.3f\n', opt.cfg.yparam, opt.yparam(valy)));
  
  % update data in mesh
  set(opt.hs, 'FaceVertexCData',     squeeze(opt.dat(:,valx,valy)));
  set(opt.hs, 'FaceVertexAlphaData', mask);

  set(opt.vline, 'xdata', [1 1]*opt.xparam(valx));
end

% update ERF-plot
set(opt.hy,    'ylim',   opt.cfg.zlim);
set(opt.vline, 'ydata',  opt.cfg.zlim);
if ~(numel(previous_vindx)==numel(opt.vindx) && all(previous_vindx==opt.vindx))
  tmp = mean(opt.dat(opt.vindx,:,valy),1);
  set(opt.tline, 'ydata', tmp);
  %set(opt.hy,    'ylim',  [min(tmp(:)) max(tmp(:))]);
  %set(opt.vline, 'ydata', [min(tmp(:)) max(tmp(:))]);
   
  set(opt.hy,    'yaxislocation', 'right');
  set(opt.stringz, 'string', sprintf('position = [%2.1f, %2.1f, %2.1f]', opt.pos(opt.vindx,:)));
end

if opt.record
  opt.frame = opt.frame + 1;
  opt.movie(opt.frame) = getframe(opt.h);
end
setappdata(h, 'opt', opt);

previous_valx = valx;
previous_valy = valy;
previous_vindx = opt.vindx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_playbutton(h, eventdata)

opt = getappdata(h, 'opt');
if strcmp(get(opt.playbutton, 'string'), 'pause')
  stop(opt.t);
  set(opt.playbutton, 'string', 'play');
else
  start(opt.t);
  set(opt.playbutton, 'string', 'pause');
end
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quitbutton(h, eventdata)

opt = getappdata(h, 'opt');
opt.cleanup = 1;
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_recordbutton(h, eventdata)

opt = getappdata(h, 'opt');
if strcmp(get(opt.recordbutton, 'string'), 'stop')
  opt.record = 0;
  set(opt.recordbutton, 'string', 'record');
else
  opt.record = 1;
  set(opt.recordbutton, 'string', 'stop');
end
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_timer(obj, info, h)

opt   = getappdata(h, 'opt');
delta = opt.speed/size(opt.dat,2);
val = get(opt.sliderx, 'value');
val = val + delta;
if val>1
  val = val-1;
end
set(opt.sliderx, 'value', val);
setappdata(h, 'opt', opt);
cb_slider(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_alim(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
switch get(h, 'String')
  case '+'
    alim(alim*sqrt(2));
  case '-'
    alim(alim/sqrt(2));
end % switch
guidata(h, opt);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_getposition(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');
if strcmp(get(get(h, 'currentaxes'), 'tag'), 'timecourse')
  % get the current point
%   pos = get(opt.hy, 'currentpoint');
%   set(opt.sliderx, 'value') = pos(1);
%   
elseif strcmp(get(get(h, 'currentaxes'), 'tag'), 'mesh')
  % get the current point
  pos = get(opt.hx, 'currentpoint');
  
  % get the view angle
  v   = get(opt.hx, 'view');
  [pos2(1),pos2(2),pos2(3)] = sph2cart(v(1),v(2),1);
  
  % get the point with the smallest positive angle to the view vector
  phi    = pos2*pos';
  [~,ix] = max(phi);
  pos    = pos(ix,:);
  
  dpos      = opt.pos - pos(ones(size(opt.pos,1),1),:);
  opt.vindx = nearest(sum(dpos.^2,2),0);
  
end
setappdata(h, 'opt', opt);
cb_slider(h);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parseKeyboardEvent(eventdata);
end
% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
    set(h, 'enable', 'off');
    drawnow;
    set(h, 'enable', 'on');
end
  
h = getparent(h);
opt = getappdata(h, 'opt');

switch key
  case 'leftarrow' % change colorlim
    opt.cfg.zlim(2) = (opt.cfg.zlim(2)-opt.cfg.zlim(1))./sqrt(2);
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.zlim);
  case 'rightarrow' % change colorlim
    opt.cfg.zlim(2) = max(opt.cfg.zlim(2).*sqrt(2), opt.cfg.zlim(1));
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.zlim);
  
  case 'shift+leftarrow'
    
  case 'shift+rightarrow'
  
  case 'uparrow' % enhance threshold
    opt.threshold = opt.threshold+0.01.*max(opt.dat(:));
    setappdata(h, 'opt', opt);
  case 'downarrow' % lower threshold
    opt.threshold = opt.threshold-0.01.*max(opt.dat(:));
    setappdata(h, 'opt', opt);
  case 'shift+uparrow' % change speed
    opt.speed = opt.speed*sqrt(2);
    setappdata(h, 'opt', opt);
  case 'shift+downarrow'
    opt.speed = opt.speed/sqrt(2);
    opt.speed = max(opt.speed, 1); % should not be smaller than 1
    setappdata(h, 'opt', opt);
  case 'ctrl+uparrow' % change channel
  case 'C' % update camera position
    camlight(opt.cam(1), 'left');
    camlight(opt.cam(2), 'right');
  case 'p'
    cb_playbutton(h);
  case 'q'
    cb_quitbutton(h);
  case 'r'
    cb_recordbutton(h);
  case 's'
    % select the speed
    response = inputdlg('speed', 'specify', 1, {num2str(opt.speed)});
    if ~isempty(response)
      opt.speed = str2double(response);
      setappdata(h, 'opt', opt);
    end
  case 't'
    % select the threshold
    response = inputdlg('threshold', 'specify', 1, {num2str(opt.threshold)});
    if ~isempty(response)
      opt.threshold = str2double(response);
      setappdata(h, 'opt', opt);
    end
  case 'z'
    % select the threshold
    response = inputdlg('colorlim', 'specify', 1, {num2str(opt.cfg.zlim(2))});
    if ~isempty(response)
      opt.cfg.zlim(2) = str2double(response);
      setappdata(h, 'opt', opt);
    end
  case 'control+control'
    % do nothing
  case 'shift+shift'
    % do nothing
  case 'alt+alt'
    % do nothing
  otherwise
    setappdata(h, 'opt', opt);
    cb_help(h);
end
cb_slider(h);
uiresume(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = parseKeyboardEvent(eventdata)

key = eventdata.Key;

% handle possible numpad events (different for Windows and UNIX systems)
% NOTE: shift+numpad number does not work on UNIX, since the shift
% modifier is always sent for numpad events
if isunix()
  shiftInd = match_str(eventdata.Modifier, 'shift');
  if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
    % now we now it was a numpad keystroke (numeric character sent AND
    % shift modifier present)
    key = eventdata.Character;
    eventdata.Modifier(shiftInd) = []; % strip the shift modifier
  end
elseif ispc()
  if strfind(eventdata.Key, 'numpad')
    key = eventdata.Character;
  end
end

if ~isempty(eventdata.Modifier)
  key = [eventdata.Modifier{1} '+' key];
end
