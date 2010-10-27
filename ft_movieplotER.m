function ft_movieplotER(cfg, timelock)

% FT_MOVIEPLOTER makes a movie of the topographic distribution of the
% time-locked average.
%
% Use as
%   ft_movieplotER(cfg, timelock)
% where the input data is from FT_TIMELOCKANALYSIS and the configuration
% can contain
%   cfg.samperframe  = number, samples per fram (default = 1)
%   cfg.framespersec = number, frames per second (default = 5)
%   cfg.xlim         = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.zlim         = 'maxmin', 'maxabs' or [zmin zmax] (default = 'maxmin')
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

% Defaults
if ~isfield(cfg, 'xlim'),          cfg.xlim = 'maxmin';           end
if ~isfield(cfg, 'zlim'),          cfg.zlim = 'maxmin';           end
if ~isfield(cfg, 'samperframe'),   cfg.samperframe = 1;           end
if ~isfield(cfg, 'framespersec'),  cfg.framespersec = 5;          end
if ~isfield(cfg, 'framesfile'),    cfg.framesfile = 'no';         end
if ~isfield(cfg, 'inputfile'),     cfg.inputfile = [];            end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    timelock = loadvar(cfg.inputfile, 'data');
  end
end

% Checkdata
timelock = ft_checkdata(timelock, 'datatype', 'timelock');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval',  {'zlim',  'absmax',  'maxabs'});

% Read or create the layout that will be used for plotting:
lay = ft_prepare_layout(cfg, timelock);

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
elseif strcmp(cfg.zlim,'maxabs')
  zmin = -max(max(abs(datavector(:))));
  zmax = max(max(abs(datavector(:))));
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end
cfg.zlim   = [zmin zmax];

figure;
% Collect frames
for iFrame = 1:floor(size(datavector,2)/cfg.samperframe)
  data_frame = mean(datavector(:,((iFrame-1)*cfg.samperframe)+1:iFrame*cfg.samperframe),2);
  % Draw topoplot:
  ft_plot_topo(chanX, chanY, data_frame, 'mask', lay.mask, 'outline', lay.outline, 'interpmethod', 'cubic');
  %ft_plot_text(chanX, chanY, '.');
  axis off
  % plot comment
  if length(ind_COMNT)==1
    comment = sprintf('zlim=[%.3g %.3g]\ntime=[%.4g %.4g]', zmin, zmax, timelock.time(((iFrame-1)*cfg.samperframe)+1), timelock.time(iFrame*cfg.samperframe));
    ft_plot_text(X_COMNT,Y_COMNT, comment)
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

uicontrol('parent', hmov, 'units', 'normalized', 'style', 'pushbutton', 'string', 'replay', 'userdata', 1, 'position', [0.86, 0.1, 0.12, 0.05], 'backgroundcolor', [1 1 1], 'callback', @replay_cb)
uicontrol('parent', hmov, 'units', 'normalized', 'style', 'pushbutton', 'string', 'frames/s', 'userdata', 1, 'position', [0.86, 0.18, 0.12, 0.05], 'backgroundcolor', [1 1 1], 'callback', @framespersec_cb)

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
