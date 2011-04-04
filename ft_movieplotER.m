function ft_movieplotER(cfg, timelock)

% FT_MOVIEPLOTER makes a movie of the topographic distribution of the
% time-locked average.
%
% Use as
%   ft_movieplotER(cfg, timelock)
% where the input data is from FT_TIMELOCKANALYSIS and the configuration
% can contain
%   cfg.xparam       = string, parameter over which the movie unrolls (default = 'time')
%   cfg.zparam       = string, parameter that is color coded (default = 'avg')
%   cfg.xlim         = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.zlim         = 'maxmin', 'maxabs' or [zmin zmax] (default = 'maxmin')
%   cfg.samperframe  = number, samples per fram (default = 1)
%   cfg.framespersec = number, frames per second (default = 5)
%   cfg.framesfile   = 'no', no file saved, or 'sting', filename of saved frames.mat (default = 'no');
%   cfg.layout       = specification of the layout, see below
%
% The layout defines how the channels are arranged. You can specify the
% layout in a variety of ways:
%  - you can provide a pre-computed layout structure (see prepare_layout)
%  - you can give the name of an ascii layout file with extension *.lay
%  - you can give the name of an electrode file
%  - you can give an electrode definition, i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% If you do not specify any of these and the data structure contains an
% electrode or gradiometer structure, that will be used for creating a
% layout. If you want to have more fine-grained control over the layout
% of the subplots, you should create your own layout file.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.

% Copyright (C) 2009, Ingrid Nieuwenhuis
% Copyright (C) 2011, Jan-Mathijs Schoffelen, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

ft_defaults;

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval',  {'zlim',  'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% set defaults
xlim         = ft_getopt(cfg, 'xlim',         'maxmin');
zlim         = ft_getopt(cfg, 'zlim',         'maxmin');
xparam       = ft_getopt(cfg, 'xparam',       'time');
zparam       = ft_getopt(cfg, 'zparam',       'avg');
samperframe  = ft_getopt(cfg, 'samperframe',  1);
framespersec = ft_getopt(cfg, 'framespersec', 5);
framesfile   = ft_getopt(cfg, 'framesfile',   []);
inputfile    = ft_getopt(cfg, 'inputfile',    []);
mask         = ft_getopt(cfg, 'mask',         []);
interactive  = ft_getopt(cfg, 'interactive', 'no');

dointeractive = istrue(interactive);

% load optional given inputfile as data
hasdata      = (nargin>1);
hasinputfile = ~isempty(inputfile);

if hasinputfile && hasdata
  error('cfg.inputfile should not be used in conjunction with giving input data to this function');
elseif hasinputfile
  timelock = loadvar(inputfile, 'data'); 
elseif hasdata
  % nothing to be done
end

% Checkdata
timelock = ft_checkdata(timelock, 'datatype', 'timelock');

% Read or create the layout that will be used for plotting:
layout = ft_prepare_layout(cfg);

% update the configuration
cfg.xparam = xparam;
cfg.zparam = zparam;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xparam = timelock.(xparam);
zparam = timelock.(zparam);

if length(xparam)~=size(zparam,2)
  error('inconsistent size of "%s" compared to "%s"', cfg.zparam, cfg.xparam);
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

% select the channels in the data that match with the layout:
[seldat, sellay] = match_str(timelock.label, layout.label);
if isempty(seldat)
  error('labels in timelock and labels in layout do not match'); 
end

% make a subselection of the data
xparam = xparam(xbeg:xend);
zparam = zparam(seldat, xbeg:xend);
clear xbeg xend

% get the x and y coordinates and labels of the channels in the data
chanX = layout.pos(sellay,1);
chanY = layout.pos(sellay,2);

% get the z-range
if ischar(zlim) && strcmp(zlim, 'maxmin')
  zlim    = [];
  zlim(1) = min(zparam(:));
  zlim(2) = max(zparam(:));
  % update the configuration
  cfg.zlim = zlim;
elseif ischar(zlim) && strcmp(cfg.zlim,'maxabs')
  zlim     = [];
  zlim(1)  = -max(abs(datavector(:)));
  zlim(2)  =  max(abs(datavector(:)));
  cfg.zlim = zlim;
end

h = gcf;
pos = get(gcf, 'position');
set(h, 'toolbar', 'figure');

if dointeractive,
    
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
  set(ht, 'position', [20 80 100 20]);
  set(ht, 'string', 'time = ');
  set(ht, 'horizontalalignment', 'left');

  %text(0,0,  sprintf('%s = \n', cfg.xparam));
  %t = timer;
  %set(t, 'timerfcn', {@cb_timer, h}, 'period', 0.1, 'executionmode', 'fixedSpacing');

  % collect the data and the options to be used in the figure
  opt.lay   = layout;
  opt.chanX = chanX;
  opt.chanY = chanY;
  opt.tim   = xparam;
  opt.dat   = zparam;
  opt.zlim  = zlim;
  opt.speed = 1;
  opt.cfg   = cfg;
  opt.s     = s;
  opt.p     = p;
  opt.t     = t;
  if ~isempty(mask) && ischar(mask)
    opt.mask = double(getsubfield(source, mask));
  end

  [dum, hs] = ft_plot_topo(chanX, chanY, zeros(numel(chanX),1), 'mask', layout.mask, 'outline', layout.outline, 'interpmethod', 'cubic');
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
  opt.ht  = ht;

  guidata(h, opt);

  % from now it is safe to hand over the control to the callback function
  set(s, 'Callback', @cb_slider);

  % from now it is safe to hand over the control to the callback function
  set(p, 'Callback', @cb_playbutton);

else
  
  % old implementation: use variables samperframe and framepersec
  % and optionally save movie in  a file
  [dum, hs] = ft_plot_topo(chanX, chanY, zeros(numel(chanX),1), 'mask', layout.mask, 'outline', layout.outline, 'interpmethod', 'cubic');
  caxis(cfg.zlim);
  axis off;
  
  xdata = get(hs, 'xdata');
  ydata = get(hs, 'ydata');
  nanmask = get(hs, 'cdata');
  for iFrame = 1:floor(size(zparam, 2)/samperframe)

    indx       = ((iFrame-1)*samperframe+1):iFrame*samperframe;
    datavector = mean(zparam(:, indx), 2);
    datamatrix = griddata(chanX, chanY, datavector, xdata, ydata, 'cubic');
    set(hs, 'cdata',  datamatrix + nanmask);
    
    F(iFrame) = getframe;
  end
  
  if ~isempty(framesfile)
    save(framesfile, 'F');
  end
  
  movie(F, 1, framespersec);
  
  % hmov = figure;
  % axis off
  % data_movie.F = F;
  % data_movie.cfg = cfg;
  % 
  % uicontrol('parent', hmov, 'units', 'normalized', 'style', 'pushbutton', 'string', 'replay', 'userdata', 1, 'position', [0.86, 0.1, 0.12, 0.05], 'backgroundcolor', [1 1 1], 'callback', @replay_cb)
  % uicontrol('parent', hmov, 'units', 'normalized', 'style', 'pushbutton', 'string', 'frames/s', 'userdata', 1, 'position', [0.86, 0.18, 0.12, 0.05], 'backgroundcolor', [1 1 1], 'callback', @framespersec_cb)
  % 
  % movie(F,1,cfg.framespersec)
  % guidata(hmov, data_movie); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath'); % this is helpful for debugging
cfg.version.id   = '$Id$'; % this will be auto-updated by the revision control system

% add information about the Matlab version used to the configuration
cfg.version.matlab = version(); % this is helpful for debugging

if isfield(timelock, 'cfg')
  % remember the configuration details of the input data
  cfg.previous = timelock.cfg;
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

% update timer info
set(opt.ht, 'string', sprintf('%s = %1.3f\n', opt.cfg.xparam, opt.tim(val)));

% update data, interpolate and render
datamatrix = griddata(opt.chanX, opt.chanY, opt.dat(:,val), opt.xdata, opt.ydata, 'cubic');
set(opt.hs, 'cdata',  datamatrix + opt.nanmask);

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