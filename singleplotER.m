function [cfg] = singleplotER(cfg, varargin)

% singleplotER plots the event-related fields or potentials of a single channel
% or the average over multiple channels. Multiple datasets can be overlayed.
%
% Use as:
%   sinlgeplotER(cfg, data)
%   singleplotER(cfg, data1, data2, ..., dataN)
%
% The data can be an ERP/ERF produced by TIMELOCKANALYSIS, a powerspectrum 
% produced by FREQANALYSIS or a coherencespectrum produced by FREQDESCRIPTIVES. 
% If you specify multiple datasets they must contain the same channels, etc.
%
% The configuration can have the following parameters:
% cfg.xparam        = field to be plotted on x-axis (default depends on data.dimord)
%                     'time' or 'freq' 
% cfg.zparam        = field to be plotted on y-axis (default depends on data.dimord)
%                     'avg', 'powspctrm' or 'cohspctrm' 
% cfg.maskparameter = field in the first dataset to be used for masking of data 
%                     (not possible for mean over multiple channels)
% cfg.maskstyle     = style used for masking of data, 'box' or 'thickness' (default = 'box')
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                     see CHANNELSELECTION for details
% cfg.cohrefchannel = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see TIMELOCKBASELINE or FREQBASELINE
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.fontsize      = font size of title (default = 8)
% cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'no')
%                     In a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. Multiple areas 
%                     can be selected by holding down the SHIFT key.
% cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
% cfg.linestyle     = linestyle/marker type, see options of the matlab PLOT function (default = '-')
% cfg.linewidth     = linewidth in points (default = 0.5)
% cfg.graphcolor    = color(s) used for plotting the dataset(s) (default = 'brgkywrgbkywrgbkywrgbkyw')
%
% See also:
%   singleplotTFR, multiplotER, multiplotTFR, topoplotER, topoplotTFR.

%
% This function depends on TIMELOCKBASELINE which has the following options:
% cfg.baseline, documented
% cfg.channel
% cfg.blcwindow
% cfg.previous
% cfg.version
%
% This function depends on FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype

% Copyright (C) 2003-2006, Ole Jensen
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

cla

% for backward compatibility with old data structures
for i=1:length(varargin)
  varargin{i} = checkdata(varargin{i});
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

GRAPHCOLOR = ['k' cfg.graphcolor ];

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
    varargin{i} = timelockanalysis(tmpcfg, varargin{i});
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
      tempdata           = freqdescriptives(tmpcfg, tempdata);
      varargin{i}.(cfg.zparam)  = tempdata.powspctrm;
      clear tempdata
    end
  else
    for i=1:(nargin-1)
      if isfield(varargin{i}, 'crsspctrm'), varargin{i} = rmfield(varargin{i}, 'crsspctrm'); end % on the fly computation of coherence spectrum is not supported
      varargin{i} = freqdescriptives(tmpcfg, varargin{i});
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
cfg = checkconfig(cfg, 'unused',  {'cohtargetchannel'});

for k=1:length(varargin)
  % Check for unconverted coherence spectrum data:
  if (strcmp(cfg.zparam,'cohspctrm')) && (isfield(varargin{k}, 'labelcmb'))
    % A reference channel is required:
    if ~isfield(cfg,'cohrefchannel'),
      error('no reference channel specified');
    end

    if strcmp(cfg.cohrefchannel, 'gui')
      % Open a single figure with the channel layout, the user can click on a reference channel
      h = clf;
      lay = prepare_layout(cfg, varargin{1});
      cfg.layout = lay;
      plot_lay(cfg.layout, 'box', false);
      title('Select the reference channel by clicking on it...');
      % add the channel information to the figure
      info       = guidata(h);
      info.x     = lay.pos(:,1);
      info.y     = lay.pos(:,2);
      info.label = lay.label;
      guidata(h, info);
      set(gcf, 'WindowButtonUpFcn', {@select_channel, 'callback', {@select_singleplotER, cfg, varargin{:}}});
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
      varargin{k} = timelockbaseline(cfg, varargin{k});
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
  cfg = checkconfig(cfg, 'renamed', {'channelindex',  'channel'});
  cfg = checkconfig(cfg, 'renamed', {'channelname',   'channel'});
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
  if nargin > 2
    colorLabels = [colorLabels inputname(k) '=' GRAPHCOLOR(k) ' '];
  end

  % select channels
  cfg.channel = channelselection(cfg.channel, varargin{k-1}.label);
  if isempty(cfg.channel)
    error('no channels selected');
  else
    chansel = match_str(varargin{k-1}.label, cfg.channel);
  end

  % Average across selected channels:
  P = squeeze(mean(P(chansel,:), 1));
  
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
    ymin = min([ymin P(ind_xmin:ind_xmax)]);
    ymax = max([ymax P(ind_xmin:ind_xmax)]);
  else
    ymin = cfg.ylim(1);
    ymax = cfg.ylim(2);
  end

  color = GRAPHCOLOR(k);
  plot_vector(varargin{k-1}.(cfg.xparam), P, 'style', cfg.linestyle, 'color', color, 'highlight', M, 'highlightstyle', cfg.maskstyle, 'linewidth', cfg.linewidth);  
end

% Set xlim and ylim:
xlim([xmin xmax]);
ylim([ymin ymax]);

% Make the figure interactive
if strcmp(cfg.interactive, 'yes')
  set(gcf, 'WindowButtonUpFcn',     {@select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'WindowButtonMotionFcn'});
end

% Create title text containing channel name(s) and channel number(s):
if length(chansel) == 1
  t = [char(cfg.channel) ' / ' num2str(chansel) ];
else
  t = sprintf('mean(%0s)', join(',', cfg.channel));
end
h = title(t,'fontsize', cfg.fontsize);

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

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
% SUBFUNCTION which is called by select_channel in case cfg.cohrefchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label, cfg, varargin)
cfg.cohrefchannel = label;
fprintf('selected cfg.cohrefchannel = ''%s''\n', cfg.cohrefchannel);
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
singleplotER(cfg, varargin{:});

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
topoplotER(cfg, varargin{:});
