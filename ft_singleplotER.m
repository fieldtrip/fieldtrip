function [cfg] = ft_singleplotER(cfg, varargin)

% ft_singleplotER plots the event-related fields or potentials of a single channel
% or the average over multiple channels. Multiple datasets can be overlayed.
%
% Use as:
%   ft_sinlgeplotER(cfg, data)
%   ft_singleplotER(cfg, data1, data2, ..., dataN)
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
% cfg.maskparameter = field in the first dataset to be used for masking of data 
%                     (not possible for mean over multiple channels, or when input contains multiple subjects
%                     or trials)
% cfg.maskstyle     = style used for masking of data, 'box', 'thickness' or 'saturation' (default = 'box')
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
% cfg.cohrefchannel = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see FT_TIMELOCKBASELINE 
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.fontsize      = font size of title (default = 8)
% cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'no')
%                     In a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. Multiple areas 
%                     can be selected by holding down the SHIFT key.
% cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
% cfg.linestyle     = linestyle/marker type, see options of the matlab PLOT function (default = '-')
%                     can be a single style for all datasets, or a cell-array containing one style for each dataset
% cfg.linewidth     = linewidth in points (default = 0.5)
% cfg.graphcolor    = color(s) used for plotting the dataset(s) (default = 'brgkywrgbkywrgbkywrgbkyw')
%                     alternatively, colors can be specified as Nx3 matrix of RGB values
%
% See also:
%   FT_SINGLEPLOTTFR, FT_MULTIPLOTER, FT_MULTIPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR

% This function depends on FT_TIMELOCKBASELINE which has the following options:
% cfg.baseline, documented
% cfg.channel
% cfg.baselinewindow
% cfg.previous
% cfg.version

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

ft_defaults

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

cla

% for backward compatibility with old data structures
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i});
end

% set the defaults:
if ~isfield(cfg,'baseline'),      cfg.baseline = 'no';                          end
if ~isfield(cfg,'trials'),        cfg.trials = 'all';                           end
if ~isfield(cfg,'xlim'),          cfg.xlim = 'maxmin';                          end
if ~isfield(cfg,'ylim'),          cfg.ylim = 'maxmin';                          end
if ~isfield(cfg,'fontsize'),      cfg.fontsize = 8;                             end
if ~isfield(cfg,'graphcolor'),    cfg.graphcolor = 'brgkywrgbkywrgbkywrgbkyw';  end
if ~isfield(cfg,'interactive'),   cfg.interactive = 'no';                       end
if ~isfield(cfg,'renderer'),      cfg.renderer = [];                            end
if ~isfield(cfg,'maskparameter'), cfg.maskparameter = [];                       end
if ~isfield(cfg,'linestyle'),     cfg.linestyle = '-';                          end
if ~isfield(cfg,'linewidth'),     cfg.linewidth = 0.5;                          end
if ~isfield(cfg,'maskstyle'),     cfg.maskstyle = 'box';                        end

if ischar(cfg.graphcolor)
  GRAPHCOLOR = ['k' cfg.graphcolor];
elseif isnumeric(cfg.graphcolor)
  GRAPHCOLOR = [0 0 0; cfg.graphcolor];
end

% check for linestyle being a cell-array, check it's length, and lengthen it if does not have enough styles in it
if ischar(cfg.linestyle)
  cfg.linestyle = {cfg.linestyle};
end
if (nargin-1) > 1
  if (length(cfg.linestyle) < (nargin-1)) && (length(cfg.linestyle) > 1)
    error('either specify cfg.linestyle as a cell-array with one cell for each dataset, or only specify one linestyle')
  elseif (length(cfg.linestyle) < (nargin-1)) && (length(cfg.linestyle) == 1)
    tmpstyle = cfg.linestyle{1};
    cfg.linestyle = cell(nargin-1,1);
    for idataset = 1:(nargin-1)
      cfg.linestyle{idataset} = tmpstyle;
    end
  end
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

% Old style coherence plotting with cohtargetchannel is no longer supported:
cfg = ft_checkconfig(cfg, 'unused',  {'cohtargetchannel'});

% Check for unconverted coherence spectrum data or any other bivariate metric:
dimtok  = tokenize(varargin{1}.dimord, '_');
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;
for k=1:length(varargin)
  % Check for unconverted coherence spectrum data:
  if (strcmp(cfg.zparam,'cohspctrm')) && (isfield(varargin{k}, 'labelcmb')) || ...
     (isfull && isfield(varargin{k}, cfg.zparam)), 
    % A reference channel is required:
    if ~isfield(cfg,'cohrefchannel'),
      error('no reference channel specified');
    end

    if strcmp(cfg.cohrefchannel, 'gui')
      % Open a single figure with the channel layout, the user can click on a reference channel
      h = clf;
      lay = ft_prepare_layout(cfg, varargin{1});
      cfg.layout = lay;
      ft_plot_lay(cfg.layout, 'box', false);
      title('Select the reference channel by clicking on it...');
      % add the channel information to the figure
      info       = guidata(h);
      info.x     = lay.pos(:,1);
      info.y     = lay.pos(:,2);
      info.label = lay.label;
      guidata(h, info);
      set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'callback', {@select_singleplotER, cfg, varargin{:}}});
      return
    end

    if ~isfull,
      % only works explicitly with coherence FIXME
      % Convert 2-dimensional channel matrix to a single dimension:
      sel1                  = strmatch(cfg.cohrefchannel, varargin{k}.labelcmb(:,2));
      sel2                  = strmatch(cfg.cohrefchannel, varargin{k}.labelcmb(:,1));
      fprintf('selected %d channels for coherence\n', length(sel1)+length(sel2));
      varargin{k}.cohspctrm = varargin{k}.cohspctrm([sel1;sel2],:,:);
      varargin{k}.label     = [varargin{k}.labelcmb(sel1,1);varargin{k}.labelcmb(sel2,2)];
      varargin{k}.labelcmb  = varargin{k}.labelcmb([sel1;sel2],:);
      varargin{k}           = rmfield(varargin{k}, 'labelcmb');
    else
      % general solution will be dealt with below
      % FIXME don't know if all still works
    end
  end

  % Apply baseline correction:
  if ~strcmp(cfg.baseline, 'no')
    if strcmp(cfg.xparam, 'time')
      varargin{k} = ft_timelockbaseline(cfg, varargin{k});
    elseif strcmp(cfg.xparam, 'freq')
      warning('Baseline correction not possible for powerspectra');
    else 
      warning('Baseline not applied, please set cfg.xparam');
    end
  end
end

% Determine x-axis range:
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

% Pick the channel(s)
if ~isfield(cfg,'channel')
  % set the default
  cfg.channel = 'all';
  % for backward compatibility
  cfg = ft_checkconfig(cfg, 'renamed', {'channelindex',  'channel'});
  cfg = ft_checkconfig(cfg, 'renamed', {'channelname',   'channel'});
end

hold on;
colorLabels = [];
ymin = [];
ymax = [];

% Plot one line for each data set:
for k=2:nargin
  % Get data matrix:
  P = varargin{k-1}.(cfg.zparam);
  labels = getfield(varargin{k-1}, 'label');

  % User colored labels if more than one data set is plotted:
  if nargin > 2 %% FIXME: variable colorLabels is not used ??? sashae
    if ischar(GRAPHCOLOR);        colorLabels = [colorLabels inputname(k) '=' GRAPHCOLOR(k) ' '];
    elseif isnumeric(GRAPHCOLOR); colorLabels = [colorLabels inputname(k) '=' num2str(GRAPHCOLOR(k,:)) ' '];
    end
 end

  % select channels
  cfg.channel = ft_channelselection(cfg.channel, varargin{k-1}.label);
  if isempty(cfg.channel)
    error('no channels selected');
  else
    chansel = match_str(varargin{k-1}.label, cfg.channel);
  end

  % Average across selected channels:
  if length(size(P)) > 2 %chan_chan_indexing
    
    refchan = match_str(varargin{k-1}.label,cfg.cohrefchannel);
      
    if strcmp(cfg.matrixside, 'feedback')
      P = squeeze(mean(mean(P(chansel,refchan,:),2),1));
    elseif strcmp(cfg.matrixside, 'feedforward')
      P = squeeze(mean(mean(P(refchan,chansel,:),2),1));
    elseif strcmp(cfg.matrixside, 'ff-fd')
      P = squeeze(mean(mean(P(refchan,chansel,:),2),1)) - squeeze(mean(mean(P(chansel,refchan,:),2),1));
    elseif strcmp(cfg.matrixside, 'fd-ff')
      P = squeeze(mean(mean(P(chansel,refchan,:),2),1)) - squeeze(mean(mean(P(refchan,chansel,:),2),1));
    end      
  else
    P = squeeze(mean(P(chansel,:), 1));
  end
  
  % select mask
  if ~isempty(cfg.maskparameter) %&& masking
    M = varargin{1}.(cfg.maskparameter); % mask always from only first dataset
    M = squeeze(mean(M(chansel,:), 1));
  else
    M = [];
  end
  
  % Update ymin and ymax for the current data set:
  if strcmp(cfg.ylim, 'maxmin')
    ind_xmin = nearest(varargin{k-1}.(cfg.xparam), xmin);
    ind_xmax = nearest(varargin{k-1}.(cfg.xparam), xmax);
    ymin = min([ymin min(P(ind_xmin:ind_xmax))]);
    ymax = max([ymax max(P(ind_xmin:ind_xmax))]);
  else
    ymin = cfg.ylim(1);
    ymax = cfg.ylim(2);
  end

  if ischar(GRAPHCOLOR);        color = GRAPHCOLOR(k);
  elseif isnumeric(GRAPHCOLOR); color = GRAPHCOLOR(k,:);
  end
  ft_plot_vector(varargin{k-1}.(cfg.xparam), P, 'style', cfg.linestyle{k-1}, 'color', color, 'highlight', M, 'highlightstyle', cfg.maskstyle, 'linewidth', cfg.linewidth);  
end

% Set xlim and ylim:
xlim([xmin xmax]);
ylim([ymin ymax]);

% Make the figure interactive
if strcmp(cfg.interactive, 'yes')
  set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'WindowButtonMotionFcn'});
end

% Create title text containing channel name(s) and channel number(s):
if length(chansel) == 1
  t = [char(cfg.channel) ' / ' num2str(chansel) ];
else
  t = sprintf('mean(%0s)', join(',', cfg.channel));
end
h = title(t,'fontsize', cfg.fontsize);

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = join(separator,cells)
if isempty(cells)
  t = '';
  return;
end
t = char(cells{1});

for i=2:length(cells)
  t = [t separator char(cells{i})];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called by ft_select_channel in case cfg.cohrefchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label, cfg, varargin)
cfg.cohrefchannel = label;
fprintf('selected cfg.cohrefchannel = ''%s''\n', cfg.cohrefchannel);
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_singleplotER(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotER(range, cfg, varargin)
cfg.comment = 'auto';
cfg.yparam = [];
cfg.xlim = range(1:2);
fprintf('selected cfg.xlim = [%f %f]\n', cfg.xlim(1), cfg.xlim(2));
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_topoplotER(cfg, varargin{:});
