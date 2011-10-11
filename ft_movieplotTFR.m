function ft_movieplotTFR(cfg, data)

% ft_movieplottfr makes a movie of a
% time frequency representation of power or coherence that was computed
% using the ft_freqanalysis or ft_freqdescriptives functions.
%
% use as
%   ft_movieplottfr(cfg, data)
% the configuration can contain
%   cfg.parameter    = string, parameter that is color coded (default = 'avg')
%   cfg.xlim         = selection boundaries over first dimension in data (e.g., time)
%                          'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim         = selection boundaries over second dimension in data (e.g., freq)
%                          'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.zlim         = plotting limits for color dimension, 'maxmin',
%                          'maxabs' or [zmin zmax] (default = 'maxmin')
%   cfg.samperframe  = number, samples per fram (default = 1)
%   cfg.framespersec = number, frames per second (default = 5)
%   cfg.framesfile   = [] (optional), no file saved, or 'string', filename of saved frames.mat (default = []);
%   cfg.moviefreq    = number, movie frames are all time points at the fixed frequency moviefreq (default = []);
%   cfg.movietime    = number, movie frames are all frequencies at the fixed time movietime (default = []);
%   cfg.layout       = specification of the layout, see below
%
% the layout defines how the channels are arranged. you can specify the
% layout in a variety of ways:
%  - you can provide a pre-computed layout structure (see prepare_layout)
%  - you can give the name of an ascii layout file with extension *.mat
%  - you can give the name of an electrode file
%  - you can give an electrode definition, i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% if you do not specify any of these and the data structure contains an
% electrode or gradiometer structure, that will be used for creating a
% layout. if you want to have more fine-grained control over the layout
% of the subplots, you should create your own layout file.
%
% to facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% if you specify this option the input data will be read from a *.mat
% file on disk. this mat files should contain only a single variable named 'data',
% corresponding to the input structure.

% copyright (c) 2009, ingrid nieuwenhuis
% copyright (c) 2011, jan-mathijs schoffelen, robert oostenveld, cristiano micheli
%
% this file is part of fieldtrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    fieldtrip is free software: you can redistribute it and/or modify
%    it under the terms of the gnu general public license as published by
%    the free software foundation, either version 3 of the license, or
%    (at your option) any later version.
%
%    fieldtrip is distributed in the hope that it will be useful,
%    but without any warranty; without even the implied warranty of
%    merchantability or fitness for a particular purpose.  see the
%    gnu general public license for more details.
%
%    you should have received a copy of the gnu general public license
%    along with fieldtrip. if not, see <http://www.gnu.org/licenses/>.
%
% $id: ft_movieploter.m 4354 2011-10-05 15:06:02z crimic $

ft_defaults

% record start time and total processing time
ftfunctimer = tic();
ftfuncclock = clock();
ftfuncmem   = memtic();

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'renamedval',  {'zlim',  'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'renamed',	 {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'deprecated',  {'xparam'});

% set defaults
xlim          = ft_getopt(cfg, 'xlim', 'maxmin');
ylim          = ft_getopt(cfg, 'ylim', 'maxmin');
zlim          = ft_getopt(cfg, 'zlim', 'maxmin');
xparam        = ft_getopt(cfg, 'xparam','time');
yparam        = ft_getopt(cfg, 'yparam');                 % default is dealt with below
parameter     = ft_getopt(cfg, 'parameter', 'powspctrm'); % use power as default
inputfile     = ft_getopt(cfg, 'inputfile',    []);
samperframe   = ft_getopt(cfg, 'samperframe',  1);
framespersec  = ft_getopt(cfg, 'framespersec', 5);
framesfile    = ft_getopt(cfg, 'framesfile',   []);
moviefreq     = ft_getopt(cfg, 'moviefreq', []);
movietime     = ft_getopt(cfg, 'movietime', []);
movierpt      = ft_getopt(cfg, 'movierpt', 1);
interactive   = ft_getopt(cfg, 'interactive', 'yes');
dointeractive = istrue(interactive);
% yparam default
if isempty(yparam) && isfield(data, 'freq')
  % the default is freq (if present)
  yparam = 'freq';
end

% load optional given inputfile as data
hasdata      = (nargin>1);
hasinputfile = ~isempty(inputfile);

if hasinputfile && hasdata
  error('cfg.inputfile should not be used in conjunction with giving input data to this function');
elseif hasinputfile
  data = loadvar(inputfile, 'data');
elseif hasdata
  % nothing to be done
end

% read or create the layout that will be used for plotting:
layout = ft_prepare_layout(cfg);

% update the configuration
cfg.xparam = xparam;
cfg.yparam = yparam;
cfg.parameter = parameter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xparam    = data.(xparam);
parameter = data.(parameter);
if ~isempty(yparam)
  yparam = data.(yparam);
end

% check consistency of xparam and yparam
% NOTE: i set two different defaults for the 'chan_time' and the 'chan_freq_time' case
if isfield(data,'dimord')
  if strcmp(data.dimord,'chan_freq_time')
    if length(xparam)~=size(parameter,3)
      error('inconsistent size of "%s" compared to "%s"', cfg.parameter, cfg.xparam);
    end
    if length(yparam)~=size(parameter,2)
      error('inconsistent size of "%s" compared to "%s"', cfg.parameter, cfg.yparam);
    end
  elseif strcmp(data.dimord,'chan_time')
    if length(xparam)~=size(parameter,2)
      error('inconsistent size of "%s" compared to "%s"', cfg.parameter, cfg.xparam);
    end    
  else
    error('input data is incompatible')
  end
end

if ischar(xlim) && strcmp(xlim, 'maxmin')
  xlim    = [];
  xlim(1) = min(xparam);
  xlim(2) = max(xparam);
end

xbeg = nearest(xparam, xlim(1));
xend = nearest(xparam, xlim(2));

% update the configuration
cfg.xlim = xparam([xbeg xend]);

if ~isempty(yparam)
  if ischar(ylim) && strcmp(ylim, 'maxmin')
    ylim    = [];
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


% select the channels in the data that match with the layout:
[seldat, sellay] = match_str(data.label, layout.label);
if isempty(seldat)
  error('labels in data and labels in layout do not match');
end

% make a subselection of the data
xparam = xparam(xbeg:xend);
yparam = yparam(ybeg:yend);
if isnan(yparam)
  parameter = parameter(seldat, xbeg:xend);  
else
  parameter = parameter(seldat, ybeg:yend, xbeg:xend);
end
clear xbeg xend ybeg yend

% get the x and y coordinates and labels of the channels in the data
chanx = layout.pos(sellay,1);
chany = layout.pos(sellay,2);

% get the z-range
if ischar(zlim) && strcmp(zlim, 'maxmin')
  zlim    = [];
  zlim(1) = min(parameter(:));
  zlim(2) = max(parameter(:));
  % update the configuration
  cfg.zlim = zlim;
elseif ischar(zlim) && strcmp(cfg.zlim,'maxabs')
  zlim     = [];
  zlim(1)  = -max(abs(parameter(:)));
  zlim(2)  =  max(abs(parameter(:)));
  cfg.zlim = zlim;
end

h = gcf;
pos = get(gcf, 'position');
set(h, 'toolbar', 'figure');

if dointeractive
  
  % add the gui elements for changing the speed
  p = uicontrol('style', 'text');
  set(p, 'position', [20 75 50 20]);
  set(p, 'string', 'speed')
  button_slower = uicontrol('style', 'pushbutton');
  set(button_slower, 'position', [75 75 20 20]);
  set(button_slower, 'string', '-')
  set(button_slower, 'callback', @cb_speed);
  button_faster = uicontrol('style', 'pushbutton');
  set(button_faster, 'position', [100 75 20 20]);
  set(button_faster, 'string', '+')
  set(button_faster, 'callback', @cb_speed);
  
  % add the gui elements for changing the color limits
  p = uicontrol('style', 'text');
  set(p, 'position', [20 100 50 20]);
  set(p, 'string', 'zlim')
  button_slower = uicontrol('style', 'pushbutton');
  set(button_slower, 'position', [75 100 20 20]);
  set(button_slower, 'string', '-')
  set(button_slower, 'callback', @cb_zlim);
  button_faster = uicontrol('style', 'pushbutton');
  set(button_faster, 'position', [100 100 20 20]);
  set(button_faster, 'string', '+')
  set(button_faster, 'callback', @cb_zlim);
  
  sx = uicontrol('style', 'slider');
  set(sx, 'position', [20 5 pos(3)-160 20]);
  % note that "sx" is needed further down
  
  sy = uicontrol('style', 'slider');
  set(sy, 'position', [20 30 pos(3)-160 20]);
  % note that "sy" is needed further down
  
  p = uicontrol('style', 'pushbutton');
  set(p, 'position', [20 50 50 20]);
  set(p, 'string', 'play')
  % note that "p" is needed further down
  
  hx = uicontrol('style', 'text');
  set(hx, 'position', [pos(3)-140 5 120 20]);
  set(hx, 'string', sprintf('%s = ', cfg.xparam));
  set(hx, 'horizontalalignment', 'left');
  
  hy = uicontrol('style', 'text');
  set(hy, 'position', [pos(3)-140 30 120 20]);
  set(hy, 'string', sprintf('%s = ', cfg.yparam));
  set(hy, 'horizontalalignment', 'left');
  
  if ~hasyparam
    set(hy, 'visible', 'off')
    set(sy, 'visible', 'off')
  end
  
  t = timer;
  set(t, 'timerfcn', {@cb_timer, h}, 'period', 0.1, 'executionmode', 'fixedspacing');
  
  % collect the data and the options to be used in the figure
  opt.lay   = layout;
  opt.chanx = chanx;
  opt.chany = chany;
  opt.xparam  = xparam; % freq
  opt.yparam  = yparam; % time
  opt.dat   = parameter;
  opt.zlim  = zlim;
  opt.speed = 1;
  opt.cfg   = cfg;
  opt.sx    = sx; % slider freq
  opt.sy    = sy; % slider time
  opt.p     = p;
  opt.t     = t;
  if ~hasyparam
    opt.timdim = 2;
  else
    opt.timdim = 3;
  end
  [dum, hs] = ft_plot_topo(chanx, chany, zeros(numel(chanx),1), 'mask', layout.mask, 'outline', layout.outline, 'interpmethod', 'cubic');
  caxis(cfg.zlim);
  axis off;
  
  % add sum stuff at a higher level for quicker access in the callback
  % routine
  opt.xdata   = get(hs, 'xdata');
  opt.ydata   = get(hs, 'ydata');
  opt.nanmask = get(hs, 'cdata');
  
  % add the handle to the mesh
  opt.hs  = hs;
  
  % add the text-handle to the guidata
  opt.hx  = hx;
  opt.hy  = hy;
  
  guidata(h, opt);
  
  % from now it is safe to hand over the control to the callback function
  set(sx, 'callback', @cb_slider);
  set(sy, 'callback', @cb_slider);
  % from now it is safe to hand over the control to the callback function
  set(p, 'callback', @cb_playbutton);
  
else
  % non interactive mode
  [~, hs] = ft_plot_topo(chanx, chany, zeros(numel(chanx),1), 'mask', layout.mask, 'outline', layout.outline, 'interpmethod', 'cubic');
  caxis(cfg.zlim);
  axis off;

  xdata   = get(hs, 'xdata');
  ydata   = get(hs, 'ydata');
  nanmask = get(hs, 'cdata');
  
  % frequency/time selection
  if ~isempty(yparam) && any(~isnan(yparam)) 
    if ~isempty(movietime)
      indx = movietime;
      for iFrame = 1:floor(size(parameter, 2)/samperframe)
        indy = ((iFrame-1)*samperframe+1):iFrame*samperframe;
        datavector = squeeze(mean(parameter(:, indy,indx), 2));
        datamatrix = griddata(chanx, chany, datavector, xdata, ydata, 'cubic');
        set(hs, 'cdata',  datamatrix + nanmask);
        F(iFrame) = getframe;
      end 
    elseif ~isempty(moviefreq)
      indy = moviefreq;
      for iFrame = 1:floor(size(parameter, 3)/samperframe)
        indx = ((iFrame-1)*samperframe+1):iFrame*samperframe;
        datavector = squeeze(mean(parameter(:, indy,indx), 3));
        datamatrix = griddata(chanx, chany, datavector, xdata, ydata, 'cubic');
        set(hs, 'cdata',  datamatrix + nanmask);
        F(iFrame) = getframe;
      end      
    else
      error('Either moviefreq or movietime should contain a bin number')
    end
  else
    for iFrame = 1:floor(size(parameter, 2)/samperframe)
      indx = ((iFrame-1)*samperframe+1):iFrame*samperframe;
      datavector = mean(parameter(:, indx), 2);    
      datamatrix = griddata(chanx, chany, datavector, xdata, ydata, 'cubic');
      set(hs, 'cdata',  datamatrix + nanmask);
      F(iFrame) = getframe;
    end
  end
   
  % save movie
  if ~isempty(framesfile)
    save(framesfile, 'F');
  end
  % play movie
  movie(F, movierpt, framespersec);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath'); % this is helpful for debugging
cfg.version.id   = '$id: ft_movieploter.m 4354 2011-10-05 15:06:02z crimic $'; % this will be auto-updated by the revision control system

% add information about the matlab version used to the configuration
cfg.callinfo.matlab = version();

% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftfunctimer);
cfg.callinfo.procmem  = memtoc(ftfuncmem);
cfg.callinfo.calltime = ftfuncclock;
cfg.callinfo.user = getusername(); % this is helpful for debugging
fprintf('the call to "%s" took %d seconds and an estimated %d mb\n', mfilename, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

if isfield(data, 'cfg')
  % remember the configuration details of the input data
  cfg.previous = data.cfg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_slider(h, eventdata)
  opt = guidata(h);

  xdim = opt.timdim;
  valx = get(opt.sx, 'value');
  valx = round(valx*(size(opt.dat,xdim)-1))+1;
  valx = min(valx, size(opt.dat,xdim));
  valx = max(valx, 1);
  if valx>size(opt.dat,opt.timdim)
    valx = size(opt.dat,opt.timdim)-1;
  end

if length(size(opt.dat))>2
  ydim = 2;
  valy = get(opt.sy, 'value');
  valy = round(valy*(size(opt.dat,ydim)-1))+1;
  valy = min(valy, size(opt.dat,ydim));
  valy = max(valy, 1);

  set(opt.hx, 'string', sprintf('%s = %f\n', opt.cfg.xparam, opt.xparam(valx)));
  set(opt.hy, 'string', sprintf('%s = %f\n', opt.cfg.yparam, opt.yparam(valy)));

  % update data, interpolate and render
  datamatrix = griddata(opt.chanx, opt.chany, opt.dat(:,valy,valx), opt.xdata, opt.ydata, 'cubic');
else
  set(opt.hx, 'string', sprintf('%s = %f\n', opt.cfg.xparam, opt.xparam(valx)));
  % update data, interpolate and render
  datamatrix = griddata(opt.chanx, opt.chany, opt.dat(:,valx), opt.xdata, opt.ydata, 'cubic');  
end
set(opt.hs, 'cdata',  datamatrix + opt.nanmask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
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
% subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_timer(obj, info, h)
if ~ishandle(h)
  return
end
opt = guidata(h);
delta = opt.speed/size(opt.dat,opt.timdim);
val = get(opt.sx, 'value');
val = val + delta;
% to avoid the slider to go out of range when the speed is too high
if val+delta>2
  val = get(opt.sx, 'value');
end
if val>1
  val = val-1; 
end
set(opt.sx, 'value', val);
cb_slider(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_zlim(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
switch get(h, 'string')
  case '+'
    caxis(caxis*sqrt(2));
  case '-'
    caxis(caxis/sqrt(2));
end % switch
guidata(h, opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_speed(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
switch get(h, 'string')
  case '+'
    opt.speed = opt.speed*sqrt(2);
  case '-'
    opt.speed = opt.speed/sqrt(2);
%     opt.speed = max(opt.speed, 1); % should not be smaller than 1
end % switch
guidata(h, opt);
