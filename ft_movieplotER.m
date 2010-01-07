function movieplotER(cfg, timelock)

% MOVIEPLOTER makes a movie of the topographic distribution of timelock
% datatypes over time.
%
% Use as: movieplotER(cfg, timelock)
% 
% cfg options:
%   cfg.samperframe  = number, samples per fram (default = 1)
%   cfg.framespersec = number, frames per second (default = 5)
%   cfg.xlim         = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.zlim         = 'maxmin', 'absmax' or [zmin zmax] (default = 'maxmin')
%   cfg.framesfile   = 'no', no file saved, or 'sting', filename of saved frames.mat (default = 'no');
%   cfg.layout        = specification of the layout, see below
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

% Copyright (C) 2009, Ingrid Nieuwenhuis

% Checkdata
timelock = checkdata(timelock, 'datatype', 'timelock');

% Defaults
if ~isfield(cfg, 'xlim'),          cfg.xlim = 'maxmin';           end
if ~isfield(cfg, 'zlim'),          cfg.zlim = 'maxmin';           end
if ~isfield(cfg, 'samperframe'),   cfg.samperframe = 1;           end
if ~isfield(cfg, 'framespersec'),  cfg.framespersec = 5;          end
if ~isfield(cfg, 'framesfile'),    cfg.framesfile = 'no';         end

% Read or create the layout that will be used for plotting:
lay = prepare_layout(cfg, timelock);

% Select the channels in the data that match with the layout:
[seldat, sellay] = match_str(timelock.label, lay.label);
if isempty(seldat)
  error('labels in timelock and labels in layout do not match'); 
end

% Get physical min/max range of x:
if strcmp(cfg.xlim,'maxmin')
  xmin = min(timelock.time);
  xmax = max(timelock.time);
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Replace value with the index of the nearest bin
xmin = nearest(timelock.time, xmin);
xmax = nearest(timelock.time, xmax);

% Make datavector of whole time window
datavector = timelock.avg(seldat,xmin:xmax);

% Select x and y coordinates and labels of the channels in the data
chanX = cfg.layout.pos(sellay,1);
chanY = cfg.layout.pos(sellay,2);
ind_COMNT = strmatch('COMNT', lay.label);
if length(ind_COMNT)==1
  % remember the position of the comment
  X_COMNT = cfg.layout.pos(ind_COMNT, 1);
  Y_COMNT = cfg.layout.pos(ind_COMNT, 1);
end

% Same zlim for whole movie! Get physical min/max range of z:
if strcmp(cfg.zlim,'maxmin')
  zmin = min(datavector(:));
  zmax = max(datavector(:));
elseif strcmp(cfg.zlim,'absmax')
  zmin = -max(max(abs(datavector(:))));
  zmax = max(max(abs(datavector(:))));
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end
cfg.zlim   = [zmin zmax];

% Collect frames
for iFrame = 1:floor(size(datavector,2)/cfg.samperframe)
  data_frame = mean(datavector(:,((iFrame-1)*cfg.samperframe)+1:iFrame*cfg.samperframe),2);
  % Draw topoplot:
  plot_topo(chanX, chanY, data_frame, 'mask', lay.mask, 'outline', lay.outline, 'interpmethod', 'cubic');
  %plot_text(chanX, chanY, '.');
  axis off
  % plot comment
  if length(ind_COMNT)==1
    comment = sprintf('zlim=[%.3g %.3g]\ntime=[%.4g %.4g]', zmin, zmax, timelock.time(((iFrame-1)*cfg.samperframe)+1), timelock.time(iFrame*cfg.samperframe));
    plot_text(X_COMNT,Y_COMNT, comment)
  end
  F(iFrame) = getframe;
  cla
end
close

% save frames file
if strcmp(cfg.framesfile, 'no');
else
  save(cfg.framesfile, 'F')
end
hmov = figure;
axis off
data_movie.F = F;
data_movie.cfg = cfg;

uicontrol('parent', hmov, 'units', 'normalized', 'style', 'pushbutton', 'string', 'replay', 'userdata', 1, 'position', [0.91, 0.1, 0.08, 0.05], 'backgroundcolor', [1 1 1], 'callback', @replay_cb)
uicontrol('parent', hmov, 'units', 'normalized', 'style', 'pushbutton', 'string', 'framespersec', 'userdata', 1, 'position', [0.91, 0.3, 0.08, 0.05], 'backgroundcolor', [1 1 1], 'callback', @framespersec_cb)

movie(F,1,cfg.framespersec)
guidata(hmov, data_movie);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function replay_cb(hmov, eventdata)
data_movie = guidata(hmov);
if get(hmov, 'userdata')
  movie(data_movie.F,1,data_movie.cfg.framespersec)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function framespersec_cb(hmov, eventdata)
data_movie = guidata(hmov);
if get(hmov, 'userdata')
  response = inputdlg('frames per second', 'specify', 1, {num2str(data_movie.cfg.framespersec)});
    if ~isempty(response)
      data_movie.cfg.framespersec = str2double(response{1});
    end
    guidata(hmov, data_movie);
end
