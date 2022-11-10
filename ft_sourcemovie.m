function [cfg, M] = ft_sourcemovie(cfg, source, source2)

% FT_SOURCEMOVIE displays the source reconstruction on a cortical mesh
% and allows the user to scroll through time with a movie.
%
% Use as
%   ft_sourcemovie(cfg, source)
% where the input source data is obtained from FT_SOURCEANALYSIS, or a
% a parcellated source structure (i.e. contains a brainordinate field) and
% cfg is a configuration structure that should contain
%
%   cfg.funparameter    = string, functional parameter that is color coded (default = 'avg.pow')
%   cfg.maskparameter   = string, functional parameter that is used for opacity (default = [])
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SOURCEPLOT, FT_SOURCEINTERPOLATE, FT_SOURCEPARCELLATE

% Copyright (C) 2011-2015, Robert Oostenveld
% Copyright (C) 2012-2014, Jorn Horschig
% Copyright (C) 2018, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
if ft_abort
  return
end
ft_preamble debug
ft_preamble loadvar source
ft_preamble provenance source


% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
source = ft_checkdata(source, 'datatype', 'source', 'feedback', 'yes');

% get the options
xlim          = ft_getopt(cfg, 'time');
ylim          = ft_getopt(cfg, 'frequency');
zlim          = ft_getopt(cfg, 'funcolorlim');
olim          = ft_getopt(cfg, 'opacitylim');                     % don't use alim as variable name
cfg.xparam    = ft_getopt(cfg, 'xparam', 'time');                 % use time as default
cfg.yparam    = ft_getopt(cfg, 'yparam');                         % default is dealt with below
cfg.funparameter  = ft_getopt(cfg, 'funparameter', 'avg.pow');        % use power as default
cfg.funcolormap   = ft_getopt(cfg, 'funcolormap',  'jet');
cfg.maskparameter = ft_getopt(cfg, 'maskparameter');

if isempty(cfg.yparam) && isfield(source, 'freq')
  % the default is freq (if present)
  cfg.yparam = 'freq';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==2
  fun = getsubfield(source, cfg.funparameter); % might be avg.pow
elseif nargin>2 && isfield(source2, 'pos')
  % in this case two conditions of data seem to have been inputted
  fun  = getsubfield(source, cfg.funparameter); % might be avg.pow
  fun2 = getsubfield(source2, cfg.funparameter);
elseif nargin>2
  % assume the first data argument to be a parcellation, and the second a
  % parcellated structure
  
  % THIS DOES NOT EXIST ANYMORE: IF PARCELLATED DATA IS IN THE INPUT, IT
  % WILL BE UNPARCELLATED AUTOMATICALLY BY FT_CHECKDATA
  ft_error('bogus');
end
if size(source.pos)~=size(fun,1)
  ft_error('inconsistent number of vertices in the cortical mesh');
end

if ~isfield(source, 'tri')
  ft_error('source.tri missing, this function requires a triangulated cortical sheet as source model');
end

if ~isempty(cfg.maskparameter) && ischar(cfg.maskparameter)
  mask = double(getsubfield(source, cfg.maskparameter));
else
  mask = 0.5*ones(size(fun));
end

xparam = source.(cfg.xparam);
if length(xparam)~=size(fun,2)
  error('inconsistent size of "%s" compared to "%s"', cfg.funparameter, cfg.xparam);
end

if ~isempty(cfg.yparam)
  yparam = source.(cfg.yparam);
  if length(cfg.yparam)~=size(fun,3)
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

if ~isempty(cfg.yparam)
  if isempty(ylim)
    ylim(1) = min(yparam);
    ylim(2) = max(yparam);
  end
  ybeg = nearest(yparam, ylim(1));
  yend = nearest(yparam, ylim(2));
  % update the configuration
  cfg.ylim = yparam([ybeg yend]);
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
if exist('fun2', 'var')
  fun2 = fun2(:,xbeg:xend,ybeg:yend);
end
mask    = mask(:,xbeg:xend,ybeg:yend);
clear xbeg xend ybeg yend

if isempty(zlim)
  zlim(1) = min(fun(:));
  zlim(2) = max(fun(:));
  % update the configuration
  cfg.funcolorlim = zlim;
end

if isempty(olim)
  olim(1) = min(mask(:));
  olim(2) = max(mask(:));
  if olim(1)==olim(2)
    olim(1) = 0;
    olim(2) = 1;
  end
  % update the configuration
  cfg.opacitylim = olim;
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
opt.tri     = source.tri;
opt.vindx   = source.inside(:);
opt.speed   = 1;
opt.record  = 0;
opt.threshold = 0;
opt.frame   = 0;
opt.cleanup = false;
if isfield(source, 'parcellation')
  opt.parcellation = source.parcellation;
  opt.parcellationlabel = source.parcellationlabel;
end


% add functional data of optional third input to the opt structure
% FIXME here we should first check whether the meshes correspond!
if nargin>2 && isfield(source2, 'pos')
  opt.dat2 = fun2;
end

% get a handle to a figure
h  = figure;
set(h, 'color', [1 1 1]);
set(h, 'toolbar', 'figure');
set(h, 'visible', 'on');
set(h, 'CloseRequestFcn', @cb_quitbutton);
set(h, 'position', [100 50 1000 600]);
set(h, 'windowbuttondownfcn', @cb_getposition);

% get timer object
t = timer;
set(t, 'timerfcn', {@cb_timer, h}, 'period', 0.1, 'executionmode', 'fixedSpacing');

% make the user interface elements
cambutton    = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'light', 'userdata', 'C');
playbutton   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'play',   'userdata', 'p');
recordbutton = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'record', 'userdata', 'r');
quitbutton   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'quit',   'userdata', 'q');

thrminmin    = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'downarrow');
thrminplus   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+downarrow');
thr          = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'threshold', 'userdata', 't');
thrplusmin   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+uparrow');
thrplusplus  = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'uparrow');

spdmin       = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'shift+downarrow');
spd          = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'speed','userdata', 's');
spdplus      = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'shift+uparrow');
climminmin   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'leftarrow');
climminplus  = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+leftarrow');
clim         = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'colorlim', 'userdata', 'z');
climplusmin  = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+rightarrow');
climplusplus = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'rightarrow');
sliderx      = uicontrol('parent', h, 'units', 'normalized', 'style', 'slider',     'string', sprintf('%s = ', cfg.xparam));
stringx      = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');
slidery      = uicontrol('parent', h, 'units', 'normalized', 'style', 'slider',     'string', sprintf('%s = ', cfg.yparam));
stringy      = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');
stringz      = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');

set(cambutton,    'position', [0.095 0.28  0.09 0.05 ], 'callback', @cb_keyboard);
set(quitbutton,   'position', [0.005 0.28  0.09 0.05 ], 'callback', @cb_keyboard);
set(playbutton,   'position', [0.005 0.22  0.09 0.05 ], 'callback', @cb_keyboard);
set(recordbutton, 'position', [0.095 0.22  0.09 0.05 ], 'callback', @cb_keyboard);
set(thrminmin,    'position', [0.005 0.16  0.03 0.025], 'callback', @cb_keyboard);
set(thrminplus,   'position', [0.005 0.185 0.03 0.025], 'callback', @cb_keyboard);
set(thr,          'position', [0.035 0.16  0.12 0.05 ], 'callback', @cb_keyboard);
set(thrplusmin,   'position', [0.155 0.16  0.03 0.025], 'callback', @cb_keyboard);
set(thrplusplus,  'position', [0.155 0.185 0.03 0.025], 'callback', @cb_keyboard);
set(climminmin,   'position', [0.005 0.10  0.03 0.025], 'callback', @cb_keyboard);
set(climminplus,  'position', [0.005 0.125 0.03 0.025], 'callback', @cb_keyboard);
set(clim,         'position', [0.035 0.10  0.12 0.05 ], 'callback', @cb_keyboard);
set(climplusmin,  'position', [0.155 0.10  0.03 0.025], 'callback', @cb_keyboard);
set(climplusplus, 'position', [0.155 0.125 0.03 0.025], 'callback', @cb_keyboard);
set(spdmin,       'position', [0.005 0.04  0.03 0.05 ], 'callback', @cb_keyboard);
set(spd,          'position', [0.035 0.04  0.12 0.05 ], 'callback', @cb_keyboard);
set(spdplus,      'position', [0.155 0.04  0.03 0.05 ], 'callback', @cb_keyboard);
set(sliderx,      'position', [0.02  0.35  0.38 0.03 ], 'callback', @cb_slider);
set(slidery,      'position', [0.200 0.005 0.78 0.03 ], 'callback', @cb_slider);
set(stringx,      'position', [0.750 0.90  0.18 0.03 ]);
set(stringy,      'position', [0.750 0.85  0.18 0.03 ]);
set(stringz,      'position', [0.60  0.95  0.33 0.03 ]);

set(stringx, 'string', sprintf('%s = ', cfg.xparam));
set(stringy, 'string', sprintf('%s = ', cfg.yparam));
set(stringz, 'string', sprintf('location = '));
set(stringx, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);
set(stringy, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);
set(stringz, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);

% create axes object to contain the mesh
hx = axes;
set(hx, 'position', [0.4 0.08 0.6 0.8]);
set(hx, 'tag', 'mesh');
hs = ft_plot_mesh(source, 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5], 'vertexcolor', 0.*opt.dat(:,1,1), 'facealpha', 0.*opt.mask(:,1,1), 'clim', [0 1], 'alphalim', [0 1], 'alphamap', 'rampup', 'colormap', cfg.funcolormap, 'maskstyle', 'colormix');
            
lighting gouraud
material dull
cam1 = camlight('left');
cam2 = camlight('right');
caxis(cfg.funcolorlim);

% create axis object to contain a time course
hy = axes;
set(hy, 'position', [0.02 0.45 0.38 0.5]);
set(hy, 'yaxislocation', 'right');

if ~hasyparam
  tline = plot(opt.xparam, mean(opt.dat(opt.vindx,:))); hold on;
  abc = axis;
  axis([opt.xparam(1) opt.xparam(end) abc(3:4)]);
  vline = plot(opt.xparam(1)*[1 1], abc(3:4), 'r');
  hline1 = plot(opt.xparam([1 end]), [0 0], 'k');
  hline2 = plot(opt.xparam([1 end]), [0 0], 'k');
  
  
  if nargin>2 && isfield(source2, 'pos')
    tline2 = plot(opt.xparam, mean(opt.dat2(opt.vindx,:)), 'r'); hold on;
  end
  
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
opt.vline = vline; % handle to the vertical line in the ERF plot
opt.tline = tline; % handle to the ERF
if nargin>2 && isfield(source2, 'pos')
  opt.tline2 = tline2;
end
opt.hline1 = hline1; % handle for the horizontal line for upper threshold
opt.hline2 = hline2; % handle for hte horizontal line for lower threshold
opt.playbutton   = playbutton; % handle to the playbutton
opt.recordbutton = recordbutton; % handle to the recordbutton
opt.quitbutton   = quitbutton; % handle to the quitbutton
opt.threshold    = [0 0];

opt.t   = t;
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
ft_postamble debug
ft_postamble previous source
ft_postamble provenance

% add a menu to the figure, but only if the current figure does not have subplots
menu_fieldtrip(gcf, cfg, false);

if ~ft_nargout
  % don't return anything
  clear cfg
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_slider(h, eventdata)

persistent previous_valx previous_valy previous_vindx previous_clim previous_thr

if isempty(previous_valx)
  previous_valx = 0;
end
if isempty(previous_valy)
  previous_valy = 0;
end
if isempty(previous_clim)
  previous_clim = [0 1];
end
if isempty(previous_thr)
  previous_thr = [0 0];
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

% update stuff
if previous_valx~=valx || previous_valy~=valy || ~isequal(previous_clim, opt.cfg.funcolorlim) || ~isequal(previous_thr, opt.threshold)
  % update strings
  set(opt.stringx, 'string', sprintf('%s = %3.3f\n', opt.cfg.xparam, opt.xparam(valx)));
  set(opt.stringy, 'string', sprintf('%s = %3.3f\n', opt.cfg.yparam, opt.yparam(valy)));
  
  dat = opt.dat(:,valx,valy);
  mask(dat>opt.threshold(1)&dat<opt.threshold(2)) = 0;
  
  % convert the color-data + opacity into rgb for robust rendering
  bgcolor = repmat([0.5 0.5 0.5], [numel(mask) 1]);
  rgb     = bg_rgba2rgb(bgcolor, dat, opt.cfg.funcolormap, opt.cfg.funcolorlim, mask, 'rampup', opt.cfg.opacitylim);
  
  % update data in mesh
  set(opt.hs, 'FaceVertexCData', rgb, 'facecolor', 'interp');
  
  set(opt.vline, 'xdata', [1 1]*opt.xparam(valx));
end

% update ERF-plot
set(opt.hy,    'ylim',   opt.cfg.funcolorlim);
set(opt.vline, 'ydata',  opt.cfg.funcolorlim);
if ~(numel(previous_vindx)==numel(opt.vindx) && all(previous_vindx==opt.vindx))
  tmp = mean(opt.dat(opt.vindx,:,valy),1);
  set(opt.tline, 'ydata', tmp);
  %set(opt.hy,    'ylim',  [min(tmp(:)) max(tmp(:))]);
  %set(opt.vline, 'ydata', [min(tmp(:)) max(tmp(:))]);
  
  if isfield(opt, 'dat2')
    tmp = mean(opt.dat2(opt.vindx,:,valy),1);
    set(opt.tline2, 'ydata', tmp);
  end
  
  set(opt.hy,    'yaxislocation', 'right');
  if isfield(opt, 'parcellation')
    set(opt.stringz, 'string', sprintf('location = %s', opt.parcellationlabel{opt.parcellation(opt.vindx)}));
  else
    set(opt.stringz, 'string', sprintf('location = [%2.1f, %2.1f, %2.1f]', opt.pos(opt.vindx,:)));
  end
end

if opt.record
  tmp = get(opt.h, 'position');
  opt.frame = opt.frame + 1;
  opt.movie(opt.frame) = getframe(opt.h,[1 1 tmp(3:4)-1]);
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
uiresume(h);

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
  %pos = get(opt.hy, 'currentpoint');
  %set(opt.sliderx, 'value', nearest(opt.xparam, pos(1)));
elseif strcmp(get(get(h, 'currentaxes'), 'tag'), 'mesh')
  % get the current point, which is defined as the intersection through the
  % axis-box (in 3D)
  pos       = get(opt.hx, 'currentpoint');
  
  % get the intersection with the mesh
  [ipos, d] = intersect_line(opt.pos, opt.tri, pos(1,:), pos(2,:));
  [md, ix]  = min(abs(d));
  
  if ~isempty(ix)
    dpos      = opt.pos - ipos(ix*ones(size(opt.pos,1),1),:);
    opt.vindx = nearest(sum(dpos.^2,2),0);
  end
end
setappdata(h, 'opt', opt);
cb_slider(h);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)

if (isempty(eventdata) && ft_platform_supports('matlabversion', -Inf, '2014a')) || isa(eventdata, 'matlab.ui.eventdata.ActionData')
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
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
    cval = opt.cfg.funcolorlim(1);
    if cval<0
      cval = cval.*1.1;
    else
      cval = cval./1.1;
    end
    opt.cfg.funcolorlim(1) = cval;
    setappdata(h, 'opt', opt);
    cb_slider(h);
    set(opt.hx, 'Clim', opt.cfg.funcolorlim);
  case 'rightarrow' % change colorlim
    cval = opt.cfg.funcolorlim(2);
    if cval<0
      cval = cval.*1.1;
    else
      cval = cval./1.1;
    end
    opt.cfg.funcolorlim(2) = cval;
    setappdata(h, 'opt', opt);
    cb_slider(h);
    set(opt.hx, 'Clim', opt.cfg.funcolorlim);
  case 'shift+leftarrow'
    cval = opt.cfg.funcolorlim(1);
    if cval<0
      cval = cval./1.1;
    else
      cval = cval.*1.1;
    end
    opt.cfg.funcolorlim(1) = cval;
    setappdata(h, 'opt', opt);
    cb_slider(h);
    set(opt.hx, 'Clim', opt.cfg.funcolorlim);
  case 'shift+rightarrow'
    cval = opt.cfg.funcolorlim(2);
    if cval<0
      cval = cval./1.1;
    else
      cval = cval.*1.1;
    end
    opt.cfg.funcolorlim(2) = cval;
    setappdata(h, 'opt', opt);
    cb_slider(h);
    set(opt.hx, 'Clim', opt.cfg.funcolorlim);
  case 'uparrow' % enhance threshold
    thrval = opt.threshold(2);
    if thrval==0 && opt.cfg.funcolorlim(2)>0
      opt.threshold(2) = opt.cfg.funcolorlim(2).*0.1;
    elseif thrval==0
      opt.threshold(2) = max(opt.dat(:)).*0.1;
    else
      opt.threshold(2) = thrval.*1.1;
    end
    setappdata(h, 'opt', opt);
    cb_slider(h);
    set(opt.hline1, 'YData', [1 1].*opt.threshold(2));
  case 'shift+uparrow' % enhance threshold
    thrval = opt.threshold(2);
    if thrval>0
      opt.threshold(2) = thrval.*0.9;
    end
    setappdata(h, 'opt', opt);
    cb_slider(h);
    set(opt.hline1, 'YData', [1 1].*opt.threshold(2));
  case 'downarrow' % lower threshold
    thrval = opt.threshold(1);
    if thrval==0 && opt.cfg.funcolorlim(1)<0
      opt.threshold(1) = opt.cfg.funcolorlim(1).*0.1;
    elseif thrval==0
      opt.threshold(1) = min(opt.dat(:)).*0.1;
    else
      opt.threshold(1) = thrval.*1.1;
    end
    setappdata(h, 'opt', opt);
    cb_slider(h);
    set(opt.hline2, 'YData', [1 1].*opt.threshold(1));
  case 'shift+downarrow' % lower threshold
    thrval = opt.threshold(1);
    if thrval<0
      opt.threshold(1) = thrval.*0.9;
    end
    setappdata(h, 'opt', opt);
    cb_slider(h);
    set(opt.hline2, 'YData', [1 1].*opt.threshold(1));
  case 'ctrl+uparrow' % change speed
    opt.speed = opt.speed*sqrt(2);
    setappdata(h, 'opt', opt);
  case 'ctrl+downarrow'
    opt.speed = opt.speed/sqrt(2);
    opt.speed = max(opt.speed, 1); % should not be smaller than 1
    setappdata(h, 'opt', opt);
  case 'ctrl+rightarrow' % change channel
  case 'C' % update camera position1.373e-14
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
      tok = tokenize(response{1}, ' ');
      tok(cellfun(@isempty,tok)) = [];
      for k = 1:numel(tok)
        opt.threshold(1,k) = str2double(tok{k});
      end
      setappdata(h, 'opt', opt);
    end
    set(opt.hline1, 'YData', [1 1].*opt.threshold(2));
    set(opt.hline2, 'YData', [1 1].*opt.threshold(1));
  case 'z'
    % select the colorlim
    response = inputdlg('colorlim', 'specify', 1, {num2str(opt.cfg.funcolorlim)});
    if ~isempty(response)
      tok = tokenize(response{1}, ' ');
      tok(cellfun(@isempty,tok)) = [];
      for k = 1:2
        opt.cfg.funcolorlim(1,k) = str2double(tok{k});
      end
      setappdata(h, 'opt', opt);
      cb_slider(h);
    end
  case 'control+control'
    % do nothing
  case 'shift+shift'
    % do nothing
  case 'alt+alt'
    % do nothing
  otherwise
    setappdata(h, 'opt', opt);
    %cb_help(h);
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
