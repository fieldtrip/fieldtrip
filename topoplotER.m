function [cfg] = topoplotER(cfg, varargin)

% TOPOPLOTER plots the topographic distribution of 2-Dimensional datatypes as
% event-related fields (ERF), potentials (ERP), the powerspectrum or coherence spectum 
% that was computed using the TIMELOCKALYSIS, TIMELOCKGRANDAVERAGE, FREQANALYSIS or 
% FREQDESCRIPTIVES functions, as a 2-D circular view (looking down at the top of the head).
%
% Use as:
%   topoplotER(cfg, data)
%
% cfg.xparam        = first dimension in data in which a selection is made
%                     'time' or 'freq' (default depends on data.dimord)
% cfg.zparam        = field that contains the data to be plotted as color 
%                     'avg', 'powspctrm' or 'cohspctrm' (default depends on data.dimord)
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.zlim          = 'maxmin', 'absmax' or [zmin zmax] (default = 'maxmin')
% cfg.cohrefchannel = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see TIMELOCKBASELINE or FREQBASELINE
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.comment       = string 'no' 'auto' or 'xlim' (default = 'auto')
%                     'auto': date, xparam and zparam limits are printed
%                     'xlim': only xparam limits are printed
% cfg.commentpos    = string or two numbers, position of comment (default 'leftbottom')
%                     'lefttop' 'leftbottom' 'middletop' 'middlebottom' 'righttop' 'rightbottom'
%                     'title' to place comment as title
%                     'layout' to place comment as specified for COMNT in layout
%                     [x y] coordinates
% cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'no')
%                     In a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. Multiple areas 
%                     can be selected by holding down the SHIFT key.
% cfg.layout        = specification of the layout, see below
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
% TOPOPLOTER calls the function TOPOPLOT to do the actual plotting. See
% the help of that function for more configuration options.
%
% See also:
%   topoplot, topoplotTFR, singleplotER, multiplotER, prepare_layout

% Undocumented local options:
% cfg.layoutname
% The following additional cfg parameters are used when plotting 3-dimensional
% data (i.e. when topoplotTFR calls topoplotER):
% cfg.yparam          field to be plotted on y-axis
% cfg.ylim            'maxmin' or [ymin ymax]         (default = 'maxmin')

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

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% $Log: topoplotER.m,v $
% Revision 1.58  2009/07/14 13:51:32  roboos
% some small changes related to the new interactive data selection
%
% Revision 1.57  2009/07/14 13:21:09  roboos
% changed the interactive plotting: instead of using plotSelection it now uses the selection function from the new plotting module (select_range and select_channel) and uses a local subfunction to update the cfg and call the next figure
%
% Revision 1.56  2009/06/17 14:05:25  roboos
% use ischar instead of isstr
%
% Revision 1.55  2009/06/17 13:44:52  roboos
% cleaned up help
%
% Revision 1.54  2009/05/12 18:58:49  roboos
% added handling of cfg.cohrefchannel='gui' for manual/interactive selection
%
% Revision 1.53  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.52  2008/12/16 15:31:42  sashae
% plot functions can now give cfg as output
% added checkconfig to start and end of function, configtracking possible
%
% Revision 1.51  2008/12/02 16:09:04  sashae
% replaced backward compatibility code by call to checkconfig,
% renamed cfg options in code to make them consistent with documentation
%
% Revision 1.50  2008/11/28 21:43:27  sashae
% allow averaging over rpt/subj also for other fields than zparam=powspctrm (thanks to Jurrian)
% added call to checkconfig
%
% Revision 1.49  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.48  2008/09/22 14:24:33  roboos
% fixed bug that was visible in case of interactive=yes, removed the scaling of the channels to unit sphere, since that is now done in prepare_layout.
%
% Revision 1.47  2008/09/22 12:31:58  roboos
% changed handling of the layout
%
% Revision 1.46  2008/04/09 07:36:09  marvger
% regular update
%
% Revision 1.45  2008/01/29 19:43:33  sashae
% added option for trial selection; plot functions now also accept data with
% repetitions (either trials or subjects), the avg is computed and plotted
% removed some old code
%
% Revision 1.44  2007/06/05 16:13:33  ingnie
% added dimord rpt_chan_time
%
% Revision 1.43  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.42  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.41  2007/03/21 15:40:07  chrhes
% updated documentation regarding the fact that cfg.layout can also contain a
% layout structure obtained using the function prepare_layout.m
%
% Revision 1.40  2007/03/14 10:45:57  chrhes
% fixed a small bug before call to topoplot to do with removal of the field
% cfg.layout from tmpcfg
%
% Revision 1.39  2007/03/14 08:43:12  roboos
% replaced call to createlayout to prepare_layout, made some small changes to the lay structure
%
% Revision 1.38  2007/01/09 10:41:43  roboos
% Added the option cfg.renderer, default is opengl. Interactive plotting
% on linux/VNC sometimes does not work, using cfg.renderer='painters'
% seems to fix it.
%
% Revision 1.37  2006/07/20 18:50:18  ingnie
% minor change to axis scaling
%
% Revision 1.36  2006/07/17 12:39:08  ingnie
% added absmax option for cfg.zlim, updated documentation
%
% Revision 1.35  2006/06/19 11:11:37  roboos
% fixed small bug in the conversion of coherence data: first select labels for the channels, then for the channelcombinations
%
% Revision 1.34  2006/05/30 14:16:42  ingnie
% updated documentation
%
% Revision 1.33  2006/05/26 12:47:28  ingnie
% added error when labels in layout and labels in data do not match and therefore
% no data is selected to be plotted
%
% Revision 1.32  2006/05/23 16:05:21  ingnie
% updated documentation
%
% Revision 1.31  2006/05/09 12:21:23  ingnie
% use OldStyle way to call topoplot.m, added some comment options, added comment
% position, updated help
%
% Revision 1.30  2006/04/27 11:42:39  jansch
% removed taking the absolute in the case of coherence-spectra
%
% Revision 1.29  2006/04/27 09:37:04  ingnie
% changed comment options, fixed bug when style=blank, added dimord
% subj_chan_time, updated documentation
%
% Revision 1.28  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.27  2006/03/22 19:09:19  jansch
% removed overwriting of cfg.layout with chanX and chanY, caused unexplicable
% crashed if toggling more than once between singleplotTFR and topoplotER in
% the interactive mode.
%
% Revision 1.26  2006/03/17 14:30:38  denpas
% Removed local drawTopoplot subroutine. topoplotER.m now calls topoplot.m to
% plot the topoplot.
%
% Revision 1.25  2006/03/14 14:55:09  roboos
% fixed detection of datatype for timelockbaseline or ferqbaseline
% changed from DOS into UNIX format
%
% Revision 1.24  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'unused',  {'cohtargetchannel'});

cla

% Multiple data sets are not supported for topoplot:
if length(varargin)>1
  error('Multiple data sets are not supported for topoplotER/topoplotTFR.');
end

data = varargin{1};

% For backward compatibility with old data structures:
data = checkdata(data);

% Set other config defaults:
if ~isfield(cfg, 'xlim'),          cfg.xlim = 'maxmin';           end
if ~isfield(cfg, 'ylim'),          cfg.ylim = 'maxmin';           end
if ~isfield(cfg, 'zlim'),          cfg.zlim = 'maxmin';           end
if ~isfield(cfg, 'style'),         cfg.style = 'both';            end
if ~isfield(cfg, 'gridscale'),     cfg.gridscale = 67;            end
if ~isfield(cfg, 'interplimits'),  cfg.interplimits = 'head';     end
if ~isfield(cfg, 'interpolation'), cfg.interpolation = 'v4';      end
if ~isfield(cfg, 'contournum'),    cfg.contournum = 6;            end
if ~isfield(cfg, 'shading'),       cfg.shading = 'flat';          end
if ~isfield(cfg, 'comment'),       cfg.comment = 'auto';          end
if ~isfield(cfg, 'commentpos'),    cfg.commentpos = 'leftbottom'; end
if ~isfield(cfg, 'ecolor'),        cfg.ecolor = [0 0 0];          end
if ~isfield(cfg, 'emarker'),       cfg.emarker = 'o';             end
if ~isfield(cfg, 'emarkersize'),   cfg.emarkersize = 2;           end
if ~isfield(cfg, 'fontsize'),      cfg.fontsize = 8;              end
if ~isfield(cfg, 'headcolor'),     cfg.headcolor = [0 0 0];       end
if ~isfield(cfg, 'hlinewidth'),    cfg.hlinewidth = 2;            end
if ~isfield(cfg, 'baseline'),      cfg.baseline = 'no';           end   %to avoid warning in timelock/freqbaseline
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';            end
if ~isfield(cfg, 'interactive'),   cfg.interactive = 'no';        end
if ~isfield(cfg, 'renderer'),      cfg.renderer = 'opengl';       end

% Set x/y/zparam defaults according to data.dimord value:
if strcmp(data.dimord, 'chan_time')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';             end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='avg';          end
elseif strcmp(data.dimord, 'chan_freq')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='freq';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';             end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';    end
elseif strcmp(data.dimord, 'subj_chan_time') || strcmp(data.dimord, 'rpt_chan_time')
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  data = timelockanalysis(tmpcfg, data);
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';             end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='avg';          end
elseif strcmp(data.dimord, 'subj_chan_freq') || strcmp(data.dimord, 'rpt_chan_freq')
  if isfield(data, 'crsspctrm'), data = rmfield(data, 'crsspctrm'); end % on the fly computation of coherence spectrum is not supported
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  tmpcfg.jackknife = 'no';
  if isfield(cfg, 'zparam') && strcmp(cfg.zparam,'cohspctrm')
    % on the fly computation of coherence spectrum is not supported
  elseif isfield(cfg, 'zparam') && ~strcmp(cfg.zparam,'powspctrm')
    % freqdesctiptives will only work on the powspctrm field, hence a temporary copy of the data is needed
    tempdata.dimord    = data.dimord;
    tempdata.freq      = data.freq;
    tempdata.label     = data.label;
    tempdata.powspctrm = data.(cfg.zparam);
    tempdata.cfg       = data.cfg;
    tempdata           = freqdescriptives(tmpcfg, tempdata);
    data.(cfg.zparam)  = tempdata.powspctrm;
    clear tempdata
  else
    data = freqdescriptives(tmpcfg, data);
  end
  if ~isfield(cfg, 'xparam'),      cfg.xparam='freq';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';             end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';    end
elseif strcmp(data.dimord, 'chan_freq_time')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';         end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';    end
elseif strcmp(data.dimord, 'subj_chan_freq_time') || strcmp(data.dimord, 'rpt_chan_freq_time')
  if isfield(data, 'crsspctrm'), data = rmfield(data, 'crsspctrm'); end % on the fly computation of coherence spectrum is not supported
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  tmpcfg.jackknife = 'no';
  if isfield(cfg, 'zparam') && strcmp(cfg.zparam,'cohspctrm')
    % on the fly computation of coherence spectrum is not supported
  elseif isfield(cfg, 'zparam') && ~strcmp(cfg.zparam,'powspctrm')
    % freqdesctiptives will only work on the powspctrm field, hence a temporary copy of the data is needed
    tempdata.dimord    = data.dimord;
    tempdata.freq      = data.freq;
    tempdata.time      = data.time;
    tempdata.label     = data.label;
    tempdata.powspctrm = data.(cfg.zparam);
    tempdata.cfg       = data.cfg;
    tempdata           = freqdescriptives(tmpcfg, tempdata);
    data.(cfg.zparam)  = tempdata.powspctrm;
    clear tempdata
  else
    data = freqdescriptives(tmpcfg, data);
  end
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';         end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';    end
elseif strcmp(data.dimord, 'chan_comp')
  % Add a pseudo-axis with the component numbers:
  data.comp = 1:size(data.topo,2);
  % Rename the field with topographic label information:
  data.label = data.topolabel;
  if ~isfield(cfg, 'xparam'),      cfg.xparam='comp';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';             end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='topo';         end
end

% Read or create the layout that will be used for plotting:
lay = prepare_layout(cfg, data);
cfg.layout = lay;

% Create time-series of small topoplots:
if ~ischar(cfg.xlim) && length(cfg.xlim)>2
  % Switch off interactive mode:
  cfg.interactive = 'no';
  xlims = cfg.xlim;
  % Iteratively call topoplotER with different xlim values:
  for i=1:length(xlims)-1
    subplot(ceil(sqrt(length(xlims)-1)), ceil(sqrt(length(xlims)-1)), i);
    cfg.xlim = xlims(i:i+1);
    topoplotER(cfg, data);
  end
  return
end

% Check for unconverted coherence spectrum data:
if (strcmp(cfg.zparam,'cohspctrm')) && isfield(data, 'labelcmb')
  % A reference channel is required:
  if ~isfield(cfg,'cohrefchannel'),
    error('no reference channel specified');
  end

  if strcmp(cfg.cohrefchannel, 'gui')
    % Open a single figure with the channel layout, the user can click on a reference channel
    h = clf;
    plot_lay(lay, 'box', false);
    title('Select the reference channel by clicking on it...');
    % add the channel information to the figure
    info       = guidata(gcf);
    info.x     = lay.pos(:,1);
    info.y     = lay.pos(:,2);
    info.label = lay.label;
    guidata(h, info);
    set(gcf, 'WindowButtonUpFcn', {@select_channel, 'callback', {@select_topoplotER, cfg, data}});
    return
  end

  % Convert 2-dimensional channel matrix to a single dimension:
  sel1           = strmatch(cfg.cohrefchannel, data.labelcmb(:,2));
  sel2           = strmatch(cfg.cohrefchannel, data.labelcmb(:,1));
  fprintf('selected %d channels for coherence\n', length(sel1)+length(sel2));
  data.cohspctrm = data.cohspctrm([sel1;sel2],:,:);
  data.label     = [data.labelcmb(sel1,1);data.labelcmb(sel2,2)];
  data.labelcmb  = data.labelcmb([sel1;sel2],:);
  data           = rmfield(data, 'labelcmb');
end

% Apply baseline correction:
if ~strcmp(cfg.baseline, 'no')
  if strcmp(cfg.xparam, 'freq') || strcmp(cfg.yparam, 'freq')
    data = freqbaseline(cfg, data);
  else
    data = timelockbaseline(cfg, data);
  end
end

% Get physical min/max range of x:
if strcmp(cfg.xlim,'maxmin')
  xmin = min(getsubfield(data, cfg.xparam));
  xmax = max(getsubfield(data, cfg.xparam));
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Replace value with the index of the nearest bin
xmin = nearest(getsubfield(data, cfg.xparam), xmin);
xmax = nearest(getsubfield(data, cfg.xparam), xmax);

% Get physical min/max range of y:
if ~isempty(cfg.yparam)
  if strcmp(cfg.ylim,'maxmin')
    ymin = min(getsubfield(data, cfg.yparam));
    ymax = max(getsubfield(data, cfg.yparam));
  else
    ymin = cfg.ylim(1);
    ymax = cfg.ylim(2);
  end

  % Replace value with the index of the nearest bin:
  ymin = nearest(getsubfield(data, cfg.yparam), ymin);
  ymax = nearest(getsubfield(data, cfg.yparam), ymax);
end

% make dat structure with one value for each channel
dat = getsubfield(data, cfg.zparam);
if ~isempty(cfg.yparam),
  dat = dat(:, ymin:ymax, xmin:xmax);
  dat = nanmean(nanmean(dat, 2), 3);
else
  dat = dat(:, xmin:xmax);
  dat = nanmean(dat, 2);
end
dat = dat(:);

% Select the channels in the data that match with the layout:
[seldat, sellay] = match_str(data.label, cfg.layout.label);
if isempty(seldat)
  error('labels in data and labels in layout do not match'); 
end
datavector = dat(seldat);
% Select x and y coordinates and labels of the channels in the data
chanX = cfg.layout.pos(sellay,1);
chanY = cfg.layout.pos(sellay,2);
chanLabels = cfg.layout.label(sellay);

% Get physical min/max range of z:
if strcmp(cfg.zlim,'maxmin')
  zmin = min(datavector);
  zmax = max(datavector);
elseif strcmp(cfg.zlim,'absmax')
  zmin = -max(max(abs(datavector)));
  zmax = max(max(abs(datavector)));
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% specify the x and y coordinates of the comment as stated in the layout
if strcmp(cfg.commentpos,'layout') 
  cfg.commentpos = [];
  ind_COMNT = strmatch('COMNT', cfg.layout.label);
  cfg.commentpos(1) = cfg.layout.pos(ind_COMNT,1);
  cfg.commentpos(2) = cfg.layout.pos(ind_COMNT,2);
end

% make cfg.comment for topoplot.m
if strcmp(cfg.comment, 'no')
  cfg = rmfield(cfg,'comment');
elseif strcmp(cfg.comment, 'auto')
  comment = date;
  if ~isempty(cfg.xparam)
    if strcmp(cfg.xlim,'maxmin')
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, cfg.xparam, data.(cfg.xparam)(xmin), data.(cfg.xparam)(xmax));
    else
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, cfg.xparam, cfg.xlim(1), cfg.xlim(2));
    end
  end
  if ~isempty(cfg.yparam)
    if strcmp(cfg.ylim,'maxmin')
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, cfg.yparam, data.(cfg.yparam)(ymin), data.(cfg.yparam)(ymax));
    else
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, cfg.yparam, cfg.ylim(1), cfg.ylim(2));
    end
  end
  if ~isempty(cfg.zparam)
    comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, cfg.zparam, zmin, zmax);
  end
  cfg.comment = comment;
elseif strcmp(cfg.comment, 'xlim')
  if strcmp(cfg.xlim,'maxmin')
    comment = sprintf('%0s=[%.3g %.3g]', cfg.xparam, data.(cfg.xparam)(xmin), data.(cfg.xparam)(xmax));
  else
    comment = sprintf('%0s=[%.3g %.3g]', cfg.xparam, cfg.xlim(1), cfg.xlim(2));
  end
  cfg.comment = comment;
elseif ~ischar(cfg.comment)
  error('cfg.comment must be string');
end

% Draw topoplot:
topoplot(cfg,chanX,chanY,datavector,chanLabels);

% The remainder of the code is meant to make the figure interactive
hold on;

% Make the figure interactive
if strcmp(cfg.interactive, 'yes')
    % add the channel information to the figure
    info       = guidata(gcf);
    info.x     = lay.pos(:,1);
    info.y     = lay.pos(:,2);
    info.label = lay.label;
    guidata(gcf, info);

  if any(strcmp(data.dimord, {'chan_time', 'chan_freq', 'subj_chan_time', 'rpt_chan_time'}))
    set(gcf, 'WindowButtonUpFcn',     {@select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonMotionFcn'});
  elseif any(strcmp(data.dimord, {'chan_freq_time', 'subj_chan_freq_time', 'rpt_chan_freq_time', 'rpttap_chan_freq_time'}))
    set(gcf, 'WindowButtonUpFcn',     {@select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, varargin{:}}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, varargin{:}}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, varargin{:}}, 'event', 'WindowButtonMotionFcn'});
  else
    error('unsupported dimord "%" for interactive plotting', data.dimord);
  end
end

axis off;
hold off;

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.cohrefchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotER(label, cfg, varargin)
cfg.cohrefchannel = label;
fprintf('selected cfg.cohrefchannel = ''%s''\n', cfg.cohrefchannel);
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
topoplotER(cfg, varargin{:});

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
  singleplotER(cfg, varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotTFR(label, cfg, varargin)
if ~isempty(label)
  cfg.xlim = 'maxmin';
  cfg.ylim = 'maxmin';
  cfg.channel = label;
  fprintf('selected cfg.channel = {');
  for i=1:(length(cfg.channel)-1)
    fprintf('''%s'', ', cfg.channel{i});
  end
  fprintf('''%s''}\n', cfg.channel{end});
  p = get(gcf, 'Position');
  f = figure;
  set(f, 'Position', p);
  singleplotTFR(cfg, varargin{:});
end
