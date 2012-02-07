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
ft_preamble callinfo
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

h = gcf;
pos = get(gcf, 'position');
set(h, 'toolbar', 'figure');

% add the GUI elements for changing the speed
p = uicontrol('style', 'text');
set(p, 'position', [20 50 50 20]);
set(p, 'string', 'speed')
button_slower = uicontrol('style', 'pushbutton');
set(button_slower, 'position', [75 50 20 20]);
set(button_slower, 'string', '-')
set(button_slower, 'Callback', @cb_speed);
button_faster = uicontrol('style', 'pushbutton');
set(button_faster, 'position', [100 50 20 20]);
set(button_faster, 'string', '+')
set(button_faster, 'Callback', @cb_speed);

% add the GUI elements for changing the color limits
p = uicontrol('style', 'text');
set(p, 'position', [20 72 50 20]);
set(p, 'string', 'zlim')
button_slower = uicontrol('style', 'pushbutton');
set(button_slower, 'position', [75 72 20 20]);
set(button_slower, 'string', '-')
set(button_slower, 'Callback', @cb_zlim);
button_faster = uicontrol('style', 'pushbutton');
set(button_faster, 'position', [100 72 20 20]);
set(button_faster, 'string', '+')
set(button_faster, 'Callback', @cb_zlim);

% add the GUI elements for changing the opacity limits
p = uicontrol('style', 'text');
set(p, 'position', [20 94 50 20]);
set(p, 'string', 'alim')
button_slower = uicontrol('style', 'pushbutton');
set(button_slower, 'position', [75 94 20 20]);
set(button_slower, 'string', '-')
set(button_slower, 'Callback', @cb_alim);
button_faster = uicontrol('style', 'pushbutton');
set(button_faster, 'position', [100 94 20 20]);
set(button_faster, 'string', '+')
set(button_faster, 'Callback', @cb_alim);

sx = uicontrol('style', 'slider');
set(sx, 'position', [20 25 pos(3)-160 20]);
% note that "sx" is needed further down

sy = uicontrol('style', 'slider');
set(sy, 'position', [20 2 pos(3)-160 20]);
% note that "sy" is needed further down

p = uicontrol('style', 'pushbutton');
set(p, 'position', [20 116 50 20]);
set(p, 'string', 'play')
% note that "p" is needed further down

m = uicontrol('style', 'pushbutton');
set(m, 'position', [20 138 50 20]);
set(m, 'string', 'record');

q = uicontrol('style', 'pushbutton');
set(q, 'position', [20 160 50 20]);
set(q, 'string', 'quit');

hx = uicontrol('style', 'text');
set(hx, 'position', [pos(3)-140 25 120 20]);
set(hx, 'string', sprintf('%s = ', cfg.xparam));
set(hx, 'horizontalalignment', 'left');

hy = uicontrol('style', 'text');
set(hy, 'position', [pos(3)-140 2 120 20]);
set(hy, 'string', sprintf('%s = ', cfg.yparam));
set(hy, 'horizontalalignment', 'left');

if ~hasyparam
  set(hy, 'Visible', 'off')
  set(sy, 'Visible', 'off')
end

t = timer;
set(t, 'timerfcn', {@cb_timer, h}, 'period', 0.1, 'executionmode', 'fixedSpacing');

% collect the data and the options to be used in the figure
opt.cfg     = cfg;
opt.xparam  = xparam;
opt.yparam  = yparam;
opt.dat     = fun;
opt.mask    = mask;
opt.speed   = 1;
opt.record  = 0;
opt.frame   = 0;
opt.cleanup = false;

if isfield(source, 'sulc')
  vdat = source.sulc;
  vdat = vdat-min(vdat)+1;
  vdat = vdat./max(vdat);
  vdat = 0.8.*repmat(vdat,[1 3]);
  ft_plot_mesh(source, 'edgecolor', 'none', 'vertexcolor', vdat);
else
  ft_plot_mesh(source, 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5]);
end
lighting gouraud
set(gca, 'Position', [0.2 0.2 0.7 0.7]);

hs = ft_plot_mesh(source, 'edgecolor', 'none', 'vertexcolor', 0*opt.dat(:,1,1), 'facealpha', 0*opt.mask(:,1,1));
lighting gouraud
cam1 = camlight('left');
cam2 = camlight('right');

caxis(cfg.zlim);
alim(cfg.alim);


% remember the varous handles
opt.h   = h;
opt.hs  = hs;
opt.p   = p;
opt.t   = t;
opt.hx  = hx;
opt.hy  = hy;
opt.sx  = sx;
opt.sy  = sy;

% add all optional information to the figure, so that the callbacks can access and modify it
guidata(h, opt);

% from now it is safe to hand over the control to the callback function
set(sx, 'Callback', @cb_slider);
set(sy, 'Callback', @cb_slider);
set(p,  'Callback', @cb_playbutton);
set(m,  'Callback', @cb_recordbutton);
set(q,  'Callback', @cb_quit);

if nargout
  % wait until the user interface is closed
  %set(h, 'CloseRequestFcn', @cb_quit);
  
  while ishandle(h)
    uiwait(h);
    opt = guidata(h);
    if opt.cleanup
      if isfield(opt, 'movie')
        M = opt.movie;
      else
        M = [];
      end
      delete(h);
    end
  end
  
  
  % do the general cleanup and bookkeeping at the end of the function
  ft_postamble trackconfig
  ft_postamble callinfo
  ft_postamble previous source

end % if nargout



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_slider(h, eventdata)
opt = guidata(h);
valx = get(opt.sx, 'value');
valx = round(valx*(size(opt.dat,2)-1))+1;
valx = min(valx, size(opt.dat,2));
valx = max(valx, 1);

valy = get(opt.sy, 'value');
valy = round(valy*(size(opt.dat,3)-1))+1;
valy = min(valy, size(opt.dat,3));
valy = max(valy, 1);

%text(0, 0, sprintf('%s = %f\n', opt.cfg.xparam, opt.tim(val)));
set(opt.hx, 'string', sprintf('%s = %f\n', opt.cfg.xparam, opt.xparam(valx)));
set(opt.hy, 'string', sprintf('%s = %f\n', opt.cfg.yparam, opt.yparam(valy)));
set(opt.hs, 'FaceVertexCData', squeeze(opt.dat(:,valx,valy)));
set(opt.hs, 'FaceVertexAlphaData', squeeze(opt.mask(:,valx,valy)));

if opt.record
  opt.frame = opt.frame + 1;
  opt.movie(opt.frame) = getframe(opt.h);
  guidata(h, opt);
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
  case 'play'
    set(h, 'string', 'stop');
    start(opt.t);
  case 'stop'
    set(h, 'string', 'play');
    stop(opt.t);
end
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_recordbutton(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
switch get(h, 'string')
  case 'record'
    set(h, 'string', 'stop');
    opt.record = 1;
  case 'stop'
    set(h, 'string', 'record');
    opt.record = 0;
end
guidata(h, opt);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quit(h, eventdata)
opt = guidata(h);
opt.cleanup = 1;
guidata(h, opt);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_timer(obj, info, h)
if ~ishandle(h)
  return
end
opt = guidata(h);
delta = opt.speed/size(opt.dat,2);
val = get(opt.sx, 'value');
val = val + delta;
if val>1
  val = val-1;
end
set(opt.sx, 'value', val);
cb_slider(h);
uiresume;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_zlim(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
switch get(h, 'String')
  case '+'
    caxis(caxis*sqrt(2));
  case '-'
    caxis(caxis/sqrt(2));
end % switch
guidata(h, opt);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_speed(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
switch get(h, 'String')
  case '+'
    opt.speed = opt.speed*sqrt(2);
  case '-'
    opt.speed = opt.speed/sqrt(2);
    opt.speed = max(opt.speed, 1); % should not be smaller than 1
end % switch
guidata(h, opt);
uiresume;

