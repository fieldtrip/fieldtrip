function cfg = ft_sourcemovie(cfg, source)

% FT_SOURCEMOVIE displays the source reconstruction on a cortical mesh 
% and allows the user to scroll through time with a movie
%
% Use as
%  FT_SOURCEMOVIE(cfg, source)
% where indata is obtained from FT_SOURCEANALYSIS
% and cfg is a configuratioun structure that should contain
%
%  cfg.parameter    = string, parameter that is color coded (default = 'avg.pow')
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

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();
ftFuncMem   = memtic();

if isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  % the input data should be read from file
  if (nargin>1)
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    source = loadvar(cfg.inputfile, 'source');
  end
end

% ensure that the input data is valiud for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
source = ft_checkdata(source, 'datatype', 'source', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'renamed',	 {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'deprecated',  {'xparam'});

% get the options
xlim    = ft_getopt(cfg, 'xlim');
zlim    = ft_getopt(cfg, 'zlim');
xparam  = ft_getopt(cfg, 'xparam', 'time');     % use time as default
parameter  = ft_getopt(cfg, 'parameter', 'avg.pow');  % use power as default
mask    = ft_getopt(cfg, 'mask',   []);

% update the configuration
cfg.xparam = xparam;
cfg.parameter = parameter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xparam = source.(xparam);
parameter = getsubfield(source, parameter); % might be avg.pow

if length(xparam)~=size(parameter,2)
  error('inconsistent size of "%s" compared to "%s"', cfg.parameter, cfg.xparam);
end

if size(source.pos)~=size(parameter,1)
  error('inconsistent number of vertices in the cortical mesh', cfg.parameter, cfg.xparam);
end

if ~isfield(source, 'tri')
  error('source.tri missing, this function requires a triangulated cortical sheet as source model');
end

if isempty(xlim)
  xlim(1) = min(xparam);
  xlim(2) = max(xparam);
end

xbeg = nearest(xparam, xlim(1));
xend = nearest(xparam, xlim(2));
% update the configuration
cfg.xlim = xparam([xbeg xend]);

% make a subselection of the data
xparam = xparam(xbeg:xend);
parameter = parameter(:,xbeg:xend);
clear xbeg xend

if isempty(zlim)
  zlim(1) = min(parameter(:));
  zlim(2) = max(parameter(:));
  % update the configuration
  cfg.zlim = zlim;
end

h = gcf;
pos = get(gcf, 'position');
set(h, 'toolbar', 'figure');

s = uicontrol('style', 'slider');
set(s, 'position', [20 20 pos(3)-40 20]);

p = uicontrol('style', 'pushbutton');
set(p, 'position', [20 50 50 20]);
set(p, 'string', 'play')

button_slower = uicontrol('style', 'pushbutton');
set(button_slower, 'position', [75 50 20 20]);
set(button_slower, 'string', '-')
set(button_slower, 'Callback', @cb_slower);

button_faster = uicontrol('style', 'pushbutton');
set(button_faster, 'position', [100 50 20 20]);
set(button_faster, 'string', '+')
set(button_faster, 'Callback', @cb_faster);

ht = uicontrol('style', 'text');
set(ht, 'position', [125 50 pos(3)-145 20]);
set(ht, 'string', 'time = ');
set(ht, 'horizontalalignment', 'left');

text(0,0,  sprintf('%s = \n', cfg.xparam));
t = timer;
set(t, 'timerfcn', {@cb_timer, h}, 'period', 0.1, 'executionmode', 'fixedSpacing');


% collect the data and the options to be used in the figure
opt.pnt = source.pos;
opt.tri = source.tri;
opt.tim = xparam;
opt.dat = parameter;
opt.speed = 1;
opt.cfg = cfg;
opt.s = s;
opt.p = p;
opt.t = t;
if ~isempty(mask) && ischar(mask)
  opt.mask = double(getsubfield(source, mask));
end

ft_plot_mesh(opt, 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5]);
lighting gouraud

hs = ft_plot_mesh(opt, 'edgecolor', 'none', 'vertexcolor', 0*opt.dat(:,1), 'facealpha', 0*opt.mask(:,1));
caxis(cfg.zlim);
lighting gouraud
camlight left
camlight right

% add the handle to the mesh
opt.hs  = hs;

% add the text-handle to the mesh
opt.ht  = ht;

guidata(h, opt);

% from now it is safe to hand over the control to the callback function
set(s, 'Callback', @cb_slider);
% from now it is safe to hand over the control to the callback function
set(p, 'Callback', @cb_playbutton);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath'); % this is helpful for debugging
cfg.version.id   = '$Id$'; % this will be auto-updated by the revision control system

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.procmem  = memtoc(ftFuncMem);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername(); % this is helpful for debugging
fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', mfilename, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

if isfield(source, 'cfg')
  % remember the configuration details of the input data
  cfg.previous = source.cfg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_slider(h, eventdata)
opt = guidata(h);
val = get(opt.s, 'value');
val = round(val*(size(opt.dat,2)-1))+1;
val = min(val, size(opt.dat,2));
val = max(val, 1);

%text(0, 0, sprintf('%s = %f\n', opt.cfg.xparam, opt.tim(val)));
set(opt.ht, 'string', sprintf('%s = %f\n', opt.cfg.xparam, opt.tim(val)));
set(opt.hs, 'FaceVertexCData', opt.dat(:,val));
if isfield(opt, 'mask')
  set(opt.hs, 'FaceVertexAlphaData', opt.mask(:,val));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_timer(obj, info, h)
if ~ishandle(h)
  return
end
opt = guidata(h);
delta = opt.speed/size(opt.dat,2);
val = get(opt.s, 'value');
val = val + delta;
if val>1
  val = val-1;
end
set(opt.s, 'value', val);
cb_slider(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_faster(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
opt.speed = opt.speed*sqrt(2);
guidata(h, opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_slower(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
opt.speed = opt.speed/sqrt(2);
opt.speed = max(opt.speed, 1); % should not be smaller than 1
guidata(h, opt);
