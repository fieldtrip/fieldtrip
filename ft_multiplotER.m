function [cfg] = ft_multiplotER(cfg, varargin)

% ft_multiplotER plots the event-related fields or potentials versus time
% or of oscillatory activity (power or coherence) versus frequency. Multiple
% datasets can be overlayed.  The plots are arranged according to their
% location specified in the layout.
%
% Use as:
%   ft_multiplotER(cfg, data)
%   ft_multiplotER(cfg, data, data2, ..., dataN)
%
% The data can be an ERP/ERF produced by FT_TIMELOCKANALYSIS, a powerspectrum
% produced by FT_FREQANALYSIS or a coherencespectrum produced by FT_FREQDESCRIPTIVES.
% If you specify multiple datasets they must contain the same channels, etc.
%
% The configuration can have the following parameters:
% cfg.xparam        = field to be plotted on x-axis (default depends on data.dimord)
%                     'time' or 'freq'
% cfg.zparam        = field to be plotted on y-axis (default depends on data.dimord)
%                     'avg', 'powspctrm' or 'cohspctrm'
% cfg.maskparameter = field in the first dataset to be used for marking significant data
% cfg.maskstyle     = style used for masking of data, 'box', 'thickness' or 'saturation' (default = 'box')
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.cohrefchannel = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see FT_TIMELOCKBASELINE or FT_FREQBASELINE
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.axes          = 'yes', 'no' (default = 'yes')
%                     Draw x- and y-axes for each graph
% cfg.box           = 'yes', 'no' (default = 'no')
%                     Draw a box around each graph
% cfg.comment       = string of text (default = date + colors)
%                     Add 'comment' to graph (according to COMNT in the layout)
% cfg.showlabels    = 'yes', 'no' (default = 'no')
% cfg.showoutline   = 'yes', 'no' (default = 'no')
% cfg.fontsize      = font size of comment and labels (if present) (default = 8)
% cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'no')
%                     In a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. Multiple areas
%                     can be selected by holding down the SHIFT key.
% cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
% cfg.linestyle     = linestyle/marker type, see options of the matlab PLOT function (default = '-')
% cfg.linewidth     = linewidth in points (default = 0.5)
% cfg.graphcolor    = color(s) used for plotting the dataset(s) (default = 'brgkywrgbkywrgbkywrgbkyw')
%                     alternatively, colors can be specified as Nx3 matrix of RGB values
%
% cfg.layout        = specify the channel layout for plotting using one of
%                     the following ways:
%
% The layout defines how the channels are arranged and what the size of each
% subplot is. You can specify the layout in a variety of ways:
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
% See also:
%   FT_MULTIPLOTTFR, FT_SINGLEPLOTER, FT_SINGLEPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR,
%   FT_PREPARE_LAYOUT

% Undocumented local options:
% cfg.layoutname
% cfg.inputfile  = one can specifiy preanalysed saved data as input
%                     The data should be provided in a cell array

%
% This function depends on FT_TIMELOCKBASELINE which has the following options:
% cfg.baseline, documented
% cfg.channel
% cfg.blcwindow
% cfg.previous
% cfg.version
%
% This function depends on FT_FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype

% Copyright (C) 2003-2006, Ole Jensen
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

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

clf

% set default for inputfile
if ~isfield(cfg, 'inputfile'),  cfg.inputfile = [];    end

hasdata = nargin>1;
if ~isempty(cfg.inputfile) % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    for i=1:numel(cfg.inputfile)
      varargin{i} = loadvar(cfg.inputfile{i}, 'data'); % read datasets from array inputfile
      data = varargin{i};
    end
  end
else
  data = varargin{1};
end

% For backward compatibility with old data structures:
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i});
end

% set the defaults:
if ~isfield(cfg,'baseline'),      cfg.baseline      = 'no';                        end
if ~isfield(cfg,'trials'),        cfg.trials        = 'all';                       end
if ~isfield(cfg,'xlim'),          cfg.xlim          = 'maxmin';                    end
if ~isfield(cfg,'ylim'),          cfg.ylim          = 'maxmin';                    end
if ~isfield(cfg,'comment'),       cfg.comment       = strcat([date '\n']);         end
if ~isfield(cfg,'axes'),          cfg.axes          = 'yes';                       end
if ~isfield(cfg,'showlabels'),    cfg.showlabels    = 'no';                        end
if ~isfield(cfg,'showoutline'),   cfg.showoutline   = 'no';                        end
if ~isfield(cfg,'box'),           cfg.box           = 'no';                        end
if ~isfield(cfg,'fontsize'),      cfg.fontsize      = 8;                           end
if ~isfield(cfg,'graphcolor'),    cfg.graphcolor    = 'brgkywrgbkywrgbkywrgbkyw';  end
if ~isfield(cfg,'interactive'),   cfg.interactive   = 'no';                        end
if ~isfield(cfg,'renderer'),      cfg.renderer      = [];                          end
if ~isfield(cfg,'maskparameter'), cfg.maskparameter = [];                          end
if ~isfield(cfg,'linestyle'),     cfg.linestyle     = '-';                         end
if ~isfield(cfg,'linewidth'),     cfg.linewidth     = 0.5;                         end
if ~isfield(cfg,'maskstyle'),     cfg.maskstyle     = 'box';                       end

if ischar(cfg.graphcolor)
  GRAPHCOLOR = ['k' cfg.graphcolor];
elseif isnumeric(cfg.graphcolor)
  GRAPHCOLOR = [0 0 0; cfg.graphcolor];
end

% Set x/y/zparam defaults according to varargin{1}.dimord value:
if strcmp(varargin{1}.dimord, 'chan_time')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='avg';                   end
elseif strcmp(varargin{1}.dimord, 'chan_freq')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
elseif strcmp(varargin{1}.dimord, 'subj_chan_time') || strcmp(varargin{1}.dimord, 'rpt_chan_time')
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  for i=1:(nargin-1)
    varargin{i} = ft_timelockanalysis(tmpcfg, varargin{i});
  end
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='avg';                   end
elseif strcmp(varargin{1}.dimord, 'subj_chan_freq') || strcmp(varargin{1}.dimord, 'rpt_chan_freq')
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  tmpcfg.jackknife = 'no';
  if isfield(cfg, 'zparam') && strcmp(cfg.zparam,'cohspctrm')
    % on the fly computation of coherence spectrum is not supported
  elseif isfield(cfg, 'zparam') && ~strcmp(cfg.zparam,'powspctrm')
    % freqdesctiptives will only work on the powspctrm field, hence a temporary copy of the data is needed
    for i=1:(nargin-1)
      tempdata.dimord    = varargin{i}.dimord;
      tempdata.freq      = varargin{i}.freq;
      tempdata.label     = varargin{i}.label;
      tempdata.powspctrm = varargin{i}.(cfg.zparam);
      tempdata.cfg       = varargin{i}.cfg;
      tempdata           = ft_freqdescriptives(tmpcfg, tempdata);
      varargin{i}.(cfg.zparam)  = tempdata.powspctrm;
      clear tempdata
    end
  else
    for i=1:(nargin-1)
      if isfield(varargin{i}, 'crsspctrm'), varargin{i} = rmfield(varargin{i}, 'crsspctrm'); end % on the fly computation of coherence spectrum is not supported
      varargin{i} = ft_freqdescriptives(tmpcfg, varargin{i});
    end
  end
  if ~isfield(cfg, 'xparam'),      cfg.xparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
end

% Make sure cfg.yparam and cfg.zparam become equivalent if only one is defined:
if (isfield(cfg, 'yparam')) && (~isfield(cfg, 'zparam'))
  cfg.zparam = cfg.yparam;
elseif (~isfield(cfg, 'yparam')) && (isfield(cfg, 'zparam'))
  cfg.yparam = cfg.zparam;
end

% Old style coherence plotting with cohtargetchannel is no longer supported
cfg = checkconfig(cfg, 'unused',  {'cohtargetchannel'});

% Read or create the layout that will be used for plotting
lay = ft_prepare_layout(cfg, varargin{1});
cfg.layout = lay;

for k=1:length(varargin)
  % Check for unconverted coherence spectrum data
  if (strcmp(cfg.zparam,'cohspctrm')) && (isfield(varargin{k}, 'labelcmb'))
    % A reference channel is required:
    if ~isfield(cfg,'cohrefchannel'),
      error('no reference channel specified');
    end
    
    if strcmp(cfg.cohrefchannel, 'gui')
      % Open a single figure with the channel layout, the user can click on a reference channel
      h = clf;
      ft_plot_lay(lay, 'box', false);
      title('Select the reference channel by clicking on it...');
      % add the channel information to the figure
      info       = guidata(h);
      info.x     = lay.pos(:,1);
      info.y     = lay.pos(:,2);
      info.label = lay.label;
      guidata(h, info);
      set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'callback', {@select_multiplotER, cfg, varargin{:}}});
      return
    end
    
    % Convert 2-dimensional channel matrix to a single dimension:
    sel1                  = strmatch(cfg.cohrefchannel, varargin{k}.labelcmb(:,2));
    sel2                  = strmatch(cfg.cohrefchannel, varargin{k}.labelcmb(:,1));
    fprintf('selected %d channels for coherence\n', length(sel1)+length(sel2));
    varargin{k}.cohspctrm = varargin{k}.cohspctrm([sel1;sel2],:,:);
    varargin{k}.label     = [varargin{k}.labelcmb(sel1,1);varargin{k}.labelcmb(sel2,2)];
    varargin{k}.labelcmb  = varargin{k}.labelcmb([sel1;sel2],:);
    varargin{k}           = rmfield(varargin{k}, 'labelcmb');
  end
  
  % Apply baseline correction:
  if ~strcmp(cfg.baseline, 'no')
    if strcmp(cfg.xparam, 'time')
      varargin{k} = ft_timelockbaseline(cfg, varargin{k});
    elseif strcmp(cfg.xparam, 'freq')
      varargin{k} = ft_freqbaseline(cfg, varargin{k});
    else
      warning('Baseline not applied, please set cfg.xparam');
    end
  end
end

% Get physical x-axis range:
if strcmp(cfg.xlim,'maxmin')
  % Find maxmin throughout all varargins:
  xmin = [];
  xmax = [];
  for i=1:length(varargin)
    xmin = min([xmin varargin{i}.(cfg.xparam)]);
    xmax = max([xmax varargin{i}.(cfg.xparam)]);
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Find corresponding x-axis bins:
xidc = find(varargin{1}.(cfg.xparam) >= xmin & varargin{1}.(cfg.xparam) <= xmax);

% Align physical x-axis range to the array bins:
xmin = varargin{1}.(cfg.xparam)(xidc(1));
xmax = varargin{1}.(cfg.xparam)(xidc(end));

% Get physical y-axis range (ylim / zparam):
if strcmp(cfg.ylim,'maxmin')
  % Find maxmin throughout all varargins:
  ymin = [];
  ymax = [];
  for i=1:length(varargin)
    % Select the channels in the data that match with the layout:
    dat = [];
    dat = getsubfield(varargin{i}, cfg.zparam);
    [seldat, sellay] = match_str(varargin{k}.label, lay.label);
    if isempty(seldat)
      error('labels in data and labels in layout do not match');
    end
    data = dat(seldat,:);
    ymin = min([ymin min(min(min(data)))]);
    ymax = max([ymax max(max(max(data)))]);
  end
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% convert the layout to Ole's style of variable names
X      = lay.pos(:,1);
Y      = lay.pos(:,2);
Width  = lay.width;
Height = lay.height;
Lbl    = lay.label;

% Create empty channel coordinates and labels arrays:
chanX(1:length(Lbl)) = NaN;
chanY(1:length(Lbl)) = NaN;
chanLabels = cell(1,length(Lbl));

hold on;
colorLabels = [];

if isfield(lay, 'outline') && strcmp(cfg.showoutline, 'yes')
  for i=1:length(lay.outline)
    if ~isempty(lay.outline{i})
      tmpX = lay.outline{i}(:,1);
      tmpY = lay.outline{i}(:,2);
      h = line(tmpX, tmpY);
      set(h, 'color', 'k');
      set(h, 'linewidth', 2);
    end
  end
end

% Plot each data set:
for k=1:length(varargin)
  P          = getsubfield(varargin{k}, cfg.zparam);
  Labels     = getfield(varargin{k}, 'label');
  
  if length(varargin) > 1
    if ischar(GRAPHCOLOR);        colorLabels = [colorLabels inputname(k+1) '=' GRAPHCOLOR(k+1) '\n'];
    elseif isnumeric(GRAPHCOLOR); colorLabels = [colorLabels inputname(k+1) '=' num2str(GRAPHCOLOR(k+1,:)) '\n'];
    end
  end
  
  if ischar(GRAPHCOLOR);        color = GRAPHCOLOR(k+1);
  elseif isnumeric(GRAPHCOLOR); color = GRAPHCOLOR(k+1,:);
  end
  
  for m=1:length(Lbl)
    l = cellstrmatch(Lbl(m),Labels);
    if ~isempty(l)
      if ~isempty(cfg.maskparameter)
        mask = varargin{1}.(cfg.maskparameter)(l,:);
      else
        mask = [];
      end
      % Plot ER:
      plotWnd(varargin{k}.(cfg.xparam),P(l,:),xidc,[xmin xmax],[ymin ymax], ...
        X(m), ...
        Y(m), ...
        Width(m), ...
        Height(m), ...
        Lbl(m), ...
        cfg, color, mask) ...
        ;
      
      % Keep ER plot coordinates (at centre of ER plot), and channel labels (will be stored in the figure's UserData struct):
      chanX(m) = X(m) + 0.5 * Width(m);
      chanY(m) = Y(m) + 0.5 * Height(m);
      chanLabels{m} = Lbl{m};
    end
  end
end

% Add the colors of the different datasets to the comment:
cfg.comment = [cfg.comment colorLabels];

% Write comment text:
l = cellstrmatch('COMNT',Lbl);
if ~isempty(l)
  ft_plot_text(X(l),Y(l),sprintf(cfg.comment),'Fontsize',cfg.fontsize);
end

% Plot scales:
l = cellstrmatch('SCALE',Lbl);
if ~isempty(l)
  plotScales([xmin xmax],[ymin ymax],X(l),Y(l),Width(1),Height(1),cfg)
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  
  % add the channel information to the figure
  info       = guidata(gcf);
  info.x     = lay.pos(:,1);
  info.y     = lay.pos(:,2);
  info.label = lay.label;
  guidata(gcf, info);
  
  set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonMotionFcn'});
end

axis tight
axis off
if strcmp(cfg.box, 'yes')
  abc = axis;
  axis(abc + [-1 +1 -1 +1]*mean(abs(abc))/10)
end
orient landscape
hold off

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end


% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotScales(xlim,ylim,xpos,ypos,width,height,cfg)
x1 =  xpos;
x2 =  xpos+width;
y1 =  ypos;
y2 =  ypos+width;
ft_plot_box([xpos xpos+width ypos ypos+height],'edgecolor','b')

if xlim(1) <=  0 && xlim(2) >= 0
  xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','b');
end

if ylim(1) <= 0 && ylim(2) >= 0
  xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','b');
end

ft_plot_text( x1,y1,num2str(xlim(1),3),'rotation',90,'HorizontalAlignment','Right','VerticalAlignment','middle','Fontsize',cfg.fontsize);
ft_plot_text( x2,y1,num2str(xlim(2),3),'rotation',90,'HorizontalAlignment','Right','VerticalAlignment','middle','Fontsize',cfg.fontsize);
ft_plot_text( x2,y1,num2str(ylim(1),3),'HorizontalAlignment','Left','VerticalAlignment','bottom','Fontsize',cfg.fontsize);
ft_plot_text( x2,y2,num2str(ylim(2),3),'HorizontalAlignment','Left','VerticalAlignment','bottom','Fontsize',cfg.fontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotWnd(x,y,xidc,xlim,ylim,xpos,ypos,width,height,label,cfg,color,mask)
set(gca,'FontSize',cfg.fontsize);

x = x(xidc);
y = y(xidc);

% Clip out of bounds y values:
y(find(y > ylim(2))) = ylim(2);
y(find(y < ylim(1))) = ylim(1);

xs = xpos+width*(x-xlim(1))/(xlim(2)-xlim(1));
ys = ypos+height*(y-ylim(1))/(ylim(2)-ylim(1));

if isempty(mask) || (~isempty(mask) && strcmp(cfg.maskstyle,'box'))
  ft_plot_vector(xs, ys, 'color', color, 'style', cfg.linestyle, 'linewidth', cfg.linewidth)
elseif ~isempty(mask) && ~strcmp(cfg.maskstyle,'box') % ft_plot_vector doesnt support boxes higher than ydata yet, so a separate option remains below
  ft_plot_vector(xs, ys, 'color', color, 'style', cfg.linestyle, 'highlight', mask, 'highlightstyle', cfg.maskstyle, 'linewidth', cfg.linewidth)
end

if strcmp(cfg.showlabels,'yes')
  ft_plot_text(xpos,ypos+1.0*height,label,'Fontsize',cfg.fontsize)
end

% Draw axes:
if strcmp(cfg.axes,'yes') || strcmp(cfg.axes, 'xy')
  % Draw y axis
  xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','k')
  % Draw x axis
  xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','k')
  
elseif strcmp(cfg.axes,'x')
  % Draw x axis
  xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','k')
  
elseif strcmp(cfg.axes,'y')
  % Draw y axis
  xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','k')
end


% Draw box around plot:
if strcmp(cfg.box,'yes')
  ft_plot_box([xpos xpos+width ypos ypos+height],'edgecolor','k')
end

% Add boxes when masktyle is box, ft_plot_vector doesnt support boxes higher than ydata yet, so this code is left here
if ~isempty(mask) && strcmp(cfg.maskstyle, 'box')
  % determine how many boxes
  foundbeg = 0;
  foundend = 0;
  beg  = [];
  eind = [];
  for i = 1:length(mask)
    if ~foundbeg  && mask(i) == 1
      beg(length(beg)+1) = i;
      foundbeg = 1;
      foundend = 0;
    elseif ~foundbeg  && mask(i) == 0
      %next
    elseif ~foundend  && mask(i) == 1
      %next
    elseif ~foundend  && mask(i) == 0
      eind(length(eind)+1) = i-1;
      foundend = 1;
      foundbeg = 0;
    end
  end
  if length(eind) == length(beg)-1
    eind(length(eind)+1) = length(mask);
  end
  numbox = length(beg);
  for i = 1:numbox
    xmaskmin = xpos+width*(x(beg(i))-xlim(1))/(xlim(2)-xlim(1));
    xmaskmax = xpos+width*(x(eind(i))-xlim(1))/(xlim(2)-xlim(1));
    %plot([xmaskmin xmaskmax xmaskmax xmaskmin xmaskmin],[ypos ypos ypos+height ypos+height ypos],'r');
    hs = patch([xmaskmin xmaskmax xmaskmax xmaskmin xmaskmin],[ypos ypos ypos+height ypos+height ypos], [.6 .6 .6]);
    set(hs, 'EdgeColor', 'none');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l = cellstrmatch(str,strlist)
l = [];
for k=1:length(strlist)
  if strcmp(char(str),char(strlist(k)))
    l = [l k];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.cohrefchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_multiplotER(label, cfg, varargin)
cfg.cohrefchannel = label;
fprintf('selected cfg.cohrefchannel = ''%s''\n', cfg.cohrefchannel);
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_multiplotER(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label, cfg, varargin)
if ~isempty(label)
  cfg.xlim = 'maxmin';
  cfg.channel = label;
  fprintf('selected cfg.channel = {');
  for i=1:(length(cfg.channel)-1)
    fprintf('''%s'', ', cfg.channel{i});
  end
  fprintf('''%s''}\n', cfg.channel{end});
  p = get(gcf, 'Position');
  f = figure;
  set(f, 'Position', p);
  ft_singleplotER(cfg, varargin{:});
end

