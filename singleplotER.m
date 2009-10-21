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
% cfg.maskparameter = field in the first dataset to be used for opacity masking of data 
%                     (not possible for mean over multiple channels)
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
% cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = 'opengl')
%
% See also:
%   singleplotTFR, multiplotER, multiplotTFR, topoplotER, topoplotTFR.

% Undocumented local options:
% cfg.GraphCol
% cfg.graphcolor
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
% $Log: singleplotER.m,v $
% Revision 1.41  2009/07/14 13:51:31  roboos
% some small changes related to the new interactive data selection
%
% Revision 1.40  2009/07/14 13:21:09  roboos
% changed the interactive plotting: instead of using plotSelection it now uses the selection function from the new plotting module (select_range and select_channel) and uses a local subfunction to update the cfg and call the next figure
%
% Revision 1.39  2009/07/13 13:25:52  crimic
% inserted new plotting tool functions
%
% Revision 1.38  2009/06/17 13:44:52  roboos
% cleaned up help
%
% Revision 1.37  2009/05/12 18:21:37  roboos
% added handling of cfg.cohrefchannel='gui' for manual/interactive selection
%
% Revision 1.36  2009/02/26 10:52:12  ingnie
% made 'maxmin' scaling of zparam according to xlim
%
% Revision 1.35  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.34  2008/12/16 14:59:12  sashae
% plot functions can now give cfg as output
% added checkconfig to start and end of function, configtracking possible
%
% Revision 1.33  2008/12/16 13:17:05  sashae
% replaced backward compatibility code by call to checkconfig
% removed obsolete subfunction pwrspctm2cohspctrm
%
% Revision 1.32  2008/11/28 22:08:19  sashae
% allow averaging over rpt/subj also for other fields than zparam=powspctrm (thanks to Jurrian)
%
% Revision 1.31  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.30  2008/01/29 19:43:33  sashae
% added option for trial selection; plot functions now also accept data with
% repetitions (either trials or subjects), the avg is computed and plotted
% removed some old code
%
% Revision 1.29  2008/01/03 16:03:02  ingnie
% channel selection per dataset, to allow plotting datasets with different amount
% of channels
%
% Revision 1.28  2007/12/18 15:57:51  ingnie
% Fixed baseline correction. Should depend on xparam not yparam (thanks to Ian). Also not possible for powerspectra , therefore removed changed into warning
%
% Revision 1.27  2007/11/08 12:10:01  ingnie
% fixed cfg.fontsize, wasn't used in making title before
%
% Revision 1.26  2007/06/19 14:01:54  ingnie
% added cfg.maskparameter, changed some white spaces
%
% Revision 1.25  2007/06/05 16:13:33  ingnie
% added dimord rpt_chan_time
%
% Revision 1.24  2007/04/19 10:26:11  roboos
% Added a warning to the "Apply baseline correction" sections. If a user doesn't set yparam, baselining is not applied. (thanks to Doug)
%
% Revision 1.23  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.22  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.21  2007/01/09 10:41:43  roboos
% Added the option cfg.renderer, default is opengl. Interactive plotting
% on linux/VNC sometimes does not work, using cfg.renderer='painters'
% seems to fix it.
%
% Revision 1.20  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.19  2006/06/19 11:11:37  roboos
% fixed small bug in the conversion of coherence data: first select labels for the channels, then for the channelcombinations
%
% Revision 1.18  2006/06/07 10:24:34  ingnie
% added error when no channels are selected
%
% Revision 1.17  2006/06/06 16:24:04  ingnie
% replaced cfg.channelindex and cfg.channelname with cfg.channel for consistency
% with other functions
%
% Revision 1.16  2006/05/30 14:18:02  ingnie
% fixed bug that appeared when plotting single channel, updated documentation
%
% Revision 1.15  2006/05/19 15:39:39  ingnie
% fixed bug that appeared when cfg.channelindex is more than one channel
%
% Revision 1.14  2006/05/03 08:12:51  ingnie
% updated documentation
%
% Revision 1.13  2006/04/27 16:01:30  ingnie
% added dimord subj_chan_time
%
% Revision 1.12  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.11  2006/04/07 23:12:51  chrhes
% changed "clf" to "cla" at the beginning to allow use in conjunction with subplot in matlab
%
% Revision 1.10  2006/03/17 14:47:15  denpas
% Updated documentation.
%
% Revision 1.9  2006/03/17 14:43:10  denpas
% Fixed cfg.yparam / cfg.zparam bug. Either param may now be used to specify
% the y-axis in case of 2D data. Also implemented data.dimord awareness.
%
% Revision 1.8  2006/03/10 09:24:16  jansch
% made a fix in assigning a default to zparam, this has to be cleaned up!
%
% Revision 1.7  2006/03/02 13:54:59  jansch
% fixed multiple small bugs
%
% Revision 1.6  2006/02/28 12:43:15  roboos
% made plotting of coherence consistent between all xxxplotER functions
% made baselining consistent, use cfg.xparam to decide between freqbaseline and timelockbaseline
%
% Revision 1.5  2006/02/27 15:03:03  denpas
% many changes, most important is added interactive functionality
% made data selection consistent between different plot functions
% changed dimord for consistency
%
% Revision 1.3  2005/05/17 17:50:38  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.2  2005/04/17 18:51:05  olejen
% Option colorgraph added
%
% Revision 1.1  2005/04/06 13:59:05  olejen
% Plot EF for a single channel. Modified from multiplotER.m
%
% Revision 1.5  2005/02/07 17:12:00  roboos
% changed handling of layout files (using new function createlayout), now also supports automatic layout creation based on gradiometer/electrode definition in data, updated help, cleaned up indentation
%
% Revision 1.4  2005/01/27 09:31:49  roboos
% applied autoindentation on code, removed many empty lines and spaces,
% replaced layoutfile reading with read_lay, applied doudavs code to all input arguments
% implemented automatic detection of arguments to plot,
% updated help, changed input from p1, p2, p3... to varargin
%
% Revision 1.3  2004/09/24 15:54:54  roboos
% included the suggested improvements by Doug Davidson: added option cfg.cohtargetchannel
% and updated the help
%
% Revision 1.2  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details
%

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
if ~isfield(cfg,'graphcolor')     cfg.graphcolor = ['brgkywrgbkywrgbkywrgbkyw'];end
if ~isfield(cfg,'interactive'),   cfg.interactive = 'no';                       end
if ~isfield(cfg,'renderer'),      cfg.renderer = 'opengl';                      end
if ~isfield(cfg,'maskparameter'), cfg.maskparameter = [];                       end

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

% Convert the layout to Ole's style of variable names:
cfg.GraphCol = GRAPHCOLOR(1);

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

  style = GRAPHCOLOR(k);
  plot_vector(varargin{k-1}.(cfg.xparam), P, 'style', style, 'highlight', M);  
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
