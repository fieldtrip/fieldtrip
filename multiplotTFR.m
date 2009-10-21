function [cfg] = multiplotTFR(cfg, data)

% multiplotTFR plots time-frequency representations of power or coherence in a 
% topographical layout. The plots of the indivual sensors are arranged according 
% to their location specified in the layout.
%
% Use as:
%   multiplotTFR(cfg, data)
%
% The data can be a time-frequency representation of power or coherence that 
% was computed using the FREQANALYSIS or FREQDESCRIPTIVES functions.
%
% The configuration can have the following parameters:
% cfg.xparam           = field to be plotted on x-axis (default depends on data.dimord)
%                        'time'
% cfg.yparam           = field to be plotted on y-axis (default depends on data.dimord)
%                        'freq'
% cfg.zparam           = field to be represented as color (default depends on data.dimord)
%                        'powspctrm' or 'cohspctrm' 
% cfg.maskparameter    = field in the data to be used for opacity masking of data
% cfg.xlim             = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim             = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.zlim             = 'maxmin','absmax' or [zmin zmax] (default = 'maxmin')
% cfg.cohrefchannel    = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline         = 'yes','no' or [time1 time2] (default = 'no'), see FREQBASELINE
% cfg.baselinetype     = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials           = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.box              = 'yes', 'no' (default = 'no' if maskparameter given default = 'yes')
%                        Draw a box around each graph
% cfg.colorbar         = 'yes', 'no' (default = 'no')
% cfg.colormap         = any sized colormap, see COLORMAP
% cfg.comment          = string of text (default = date + zlimits)
%                        Add 'comment' to graph (according to COMNT in the layout)
% cfg.showlabels       = 'yes', 'no' (default = 'no')
% cfg.showoutline      = 'yes', 'no' (default = 'no')
% cfg.fontsize         = font size of comment and labels (if present) (default = 8)
% cfg.interactive      = Interactive plot 'yes' or 'no' (default = 'no')
%                        In a interactive plot you can select areas and produce a new
%                        interactive plot when a selected area is clicked. Multiple areas 
%                        can be selected by holding down the SHIFT key.
% cfg.masknans         = 'yes' or 'no' (default = 'yes')
% cfg.renderer         = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
% cfg.layout           = specify the channel layout for plotting using one of 
%                        the following ways:
%
% The layout defines how the channels are arranged and what the size of each
% subplot is. You can specify the layout in a variety of ways:
%  - you can provide a pre-computed layout structure (see prepare_layout)
%  - you can give the name of an ascii layout file with extension *.lay
%  - you can give the name of an electrode file
%  - you can give an electrode definition, i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% If you do not specify any of these and the data structure contains an
% electrode or gradiometer structure (common for MEG data, since the header
% of the MEG datafile contains the gradiometer information), that will be
% used for creating a layout. If you want to have more fine-grained control
% over the layout of the subplots, you should create your own layout file.
%
% See also:
%   multiplotER, singleplotER, singleplotTFR, topoplotER, topoplotTFR,
%   prepare_layout

% Undocumented local options:
% cfg.channel
% cfg.layoutname
% cfg.xparam
% cfg.zparam
%
% This function depends on FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype, documented

% Copyright (C) 2003-2006, Ole Jensen
%
% $Log: multiplotTFR.m,v $
% Revision 1.51  2009/07/14 13:52:16  roboos
% changed the interactive plotting: instead of using plotSelection it now uses the selection function from the new plotting module (select_range and select_channel) and uses a local subfunction to update the cfg and call the next figure
%
% Revision 1.50  2009/07/14 13:27:25  roboos
% consistent handling of cfg.renderer, default is to let matlab decide
%
% Revision 1.49  2009/06/17 13:44:52  roboos
% cleaned up help
%
% Revision 1.48  2009/05/12 18:48:15  roboos
% added handling of cfg.cohrefchannel='gui' for manual/interactive selection
% changed default for cfg.renderer -> let matlab decide
%
% Revision 1.47  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.46  2008/12/16 14:59:13  sashae
% plot functions can now give cfg as output
% added checkconfig to start and end of function, configtracking possible
%
% Revision 1.45  2008/12/16 13:17:06  sashae
% replaced backward compatibility code by call to checkconfig
% removed obsolete subfunction pwrspctm2cohspctrm
%
% Revision 1.44  2008/11/28 14:30:30  sashae
% averaging over rpt/subj skipped if zparam='cohspctrm'
%
% Revision 1.43  2008/11/27 14:48:27  roboos
% allow averaging over rpt or subj also for other fieldz than zparam=powspctrm (thanks to Jurrian)
%
% Revision 1.42  2008/10/29 12:40:58  roboos
% removed "axis equal"
%
% Revision 1.41  2008/10/28 14:30:51  ingnie
% linearity of axis is tested, warning given in case of nonlinear axis
%
% Revision 1.40  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.39  2008/09/22 12:53:27  roboos
% ensure equal and tight axes
%
% Revision 1.38  2008/09/22 12:44:08  roboos
% added option cfg.showoutline, default is no
%
% Revision 1.37  2008/09/22 12:30:59  roboos
% position the subplots with the box centered on the electrode position
%
% Revision 1.36  2008/01/29 19:43:33  sashae
% added option for trial selection; plot functions now also accept data with
% repetitions (either trials or subjects), the avg is computed and plotted
% removed some old code
%
% Revision 1.35  2007/10/12 12:58:59  marvger
% oops; removed accidental keyboard command
%
% Revision 1.34  2007/10/10 10:22:07  marvger
% *** empty log message ***
%
% Revision 1.33  2007/06/19 13:57:46  ingnie
% axis wider if cfg.box = 'yes'
%
% Revision 1.32  2007/06/14 12:23:48  ingnie
% added cfg.colormap option
%
% Revision 1.31  2007/06/13 12:09:59  ingnie
% made alphadata scaled between 0 and 1, partly transparent possible now
%
% Revision 1.30  2007/06/13 09:32:01  ingnie
% fixed maskparameter (all chans got mask of last chan, is fixed now)
%
% Revision 1.29  2007/06/06 10:04:29  jansch
% fixed typo on line 173 ~empty into ~isempty
%
% Revision 1.28  2007/06/05 16:14:58  ingnie
% added cfg.maskparameter and cfg.box
%
% Revision 1.27  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.26  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.25  2007/03/21 15:49:41  chrhes
% updated documentation regarding the fact that cfg.layout can also contain a
% layout structure obtained using the function prepare_layout.m
%
% Revision 1.24  2007/03/14 08:43:12  roboos
% replaced call to createlayout to prepare_layout, made some small changes to the lay structure
%
% Revision 1.23  2007/01/09 10:41:43  roboos
% Added the option cfg.renderer, default is opengl. Interactive plotting
% on linux/VNC sometimes does not work, using cfg.renderer='painters'
% seems to fix it.
%
% Revision 1.22  2006/06/20 16:25:59  ingnie
% updated documentation
%
% Revision 1.21  2006/06/19 11:11:37  roboos
% fixed small bug in the conversion of coherence data: first select labels for the channels, then for the channelcombinations
%
% Revision 1.20  2006/06/06 16:57:51  ingnie
% updated documentation
%
% Revision 1.19  2006/05/30 14:16:42  ingnie
% updated documentation
%
% Revision 1.18  2006/05/26 12:42:00  ingnie
% added error when labels in layout and labels in data do not match and therefore
% no data is selected to be plotted
%
% Revision 1.17  2006/05/03 08:12:51  ingnie
% updated documentation
%
% Revision 1.16  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.15  2006/03/23 09:48:35  jansch
% fixed bug in plotting of coherence-spectrum
%
% Revision 1.14  2006/03/06 11:47:55  denpas
% Fixed zmin/zmax bug.
%
% Revision 1.13  2006/02/27 15:03:03  denpas
% many changes, most important is added interactive functionality
% made data selection consistent between different plot functions
% changed dimord for consistency
%
% Revision 1.9  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.8  2005/04/08 17:26:54  olejen
% TFR was wrongly assigned before calling freqbaseline.m - this is now reversed
%
% Revision 1.7  2005/04/06 14:26:55  olejen
% clf before plotting
% now make use of freqbaseline.m
%
% Revision 1.6  2005/04/06 07:40:56  jansch
% included option cfg.cohrefchannel. updated help.
%
% Revision 1.5  2005/02/07 17:12:00  roboos
% changed handling of layout files (using new function createlayout), now also supports automatic layout creation based on gradiometer/electrode definition in data, updated help, cleaned up indentation
%
% Revision 1.4  2004/09/24 15:54:54  roboos
% included the suggested improvements by Doug Davidson: added option cfg.cohtargetchannel
% and updated the help
%
% Revision 1.3  2004/09/01 18:02:23  roboos
% added copyright statements, removed cfg as output argument
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

clf

% for backward compatibility with old data structures
data = checkdata(data);

% Set the defaults:
if ~isfield(cfg,'baseline'),        cfg.baseline = 'no';               end
if ~isfield(cfg,'baselinetype'),    cfg.baselinetype = 'absolute';     end
if ~isfield(cfg,'trials'),          cfg.trials = 'all';                end
if ~isfield(cfg,'xlim'),            cfg.xlim = 'maxmin';               end
if ~isfield(cfg,'ylim'),            cfg.ylim = 'maxmin';               end
if ~isfield(cfg,'zlim'),            cfg.zlim = 'maxmin';               end
if ~isfield(cfg,'colorbar'),        cfg.colorbar = 'no';               end
if ~isfield(cfg,'comment'),         cfg.comment = date;                end
if ~isfield(cfg,'showlabels'),      cfg.showlabels = 'no';             end
if ~isfield(cfg,'showoutline'),     cfg.showoutline = 'no';            end
if ~isfield(cfg,'channel'),         cfg.channel = 'all';               end
if ~isfield(cfg,'fontsize'),        cfg.fontsize = 8;                  end
if ~isfield(cfg,'interactive'),     cfg.interactive = 'no';            end
if ~isfield(cfg,'renderer'),        cfg.renderer = [];                 end % let matlab decide on default
if ~isfield(cfg,'masknans'),        cfg.masknans = 'yes';              end
if ~isfield(cfg,'maskparameter'),   cfg.maskparameter = [];            end
if ~isfield(cfg,'box')             
  if ~isempty(cfg.maskparameter)
    cfg.box = 'yes';
  else
    cfg.box = 'no';
  end
end

% Set x/y/zparam defaults according to data.dimord value:
if strcmp(data.dimord, 'chan_freq_time')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
elseif strcmp(data.dimord, 'subj_chan_freq_time') || strcmp(data.dimord, 'rpt_chan_freq_time')
  if isfield(data, 'crsspctrm'),  data = rmfield(data, 'crsspctrm');  end % on the fly computation of coherence spectrum is not supported
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
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
end

% Old style coherence plotting with cohtargetchannel is no longer supported:
cfg = checkconfig(cfg, 'unused',  {'cohtargetchannel'});

% Read or create the layout that will be used for plotting:
lay = prepare_layout(cfg, data);
cfg.layout = lay;

% Check for unconverted coherence spectrum data:
if (strcmp(cfg.zparam,'cohspctrm')),
  % A reference channel is required:
  if ~isfield(cfg,'cohrefchannel'),
    error('no reference channel specified');
  end

  if strcmp(cfg.cohrefchannel, 'gui')
    % Open a single figure with the channel layout, the user can click on a reference channel
    h = clf;
    plot_lay(lay, 'box', false);
    title('Select the reference channel by clicking on it...');
    info       = [];
    info.x     = lay.pos(:,1);
    info.y     = lay.pos(:,2);
    info.label = lay.label;
    guidata(h, info);
    set(gcf, 'WindowButtonUpFcn', {@select_channel, 'callback', {@select_multiplotTFR, cfg, data}});
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
  data = freqbaseline(cfg, data);
end

% Get physical x-axis range:
if strcmp(cfg.xlim,'maxmin')
  xmin = min(data.(cfg.xparam));
  xmax = max(data.(cfg.xparam));
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Find corresponding x-axis bins:
xidc = find(data.(cfg.xparam) >= xmin & data.(cfg.xparam) <= xmax);

% Align physical x-axis range to the array bins:
xmin = data.(cfg.xparam)(xidc(1));
xmax = data.(cfg.xparam)(xidc(end));

% Get physical y-axis range:
if strcmp(cfg.ylim,'maxmin')
  ymin = min(data.(cfg.yparam));
  ymax = max(data.(cfg.yparam));
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Find corresponding y-axis bins:
yidc = find(data.(cfg.yparam) >= ymin & data.(cfg.yparam) <= ymax);

% Align physical y-axis range to the array bins:
ymin = data.(cfg.yparam)(yidc(1));
ymax = data.(cfg.yparam)(yidc(end));

% test if X and Y are linearly spaced (to within 10^-12): % FROM UIMAGE
x = data.(cfg.xparam)(xidc);
y = data.(cfg.yparam)(yidc);
dx = min(diff(x));  % smallest interval for X
dy = min(diff(y));  % smallest interval for Y
evenx = all(abs(diff(x)/dx-1)<1e-12);     % true if X is linearly spaced
eveny = all(abs(diff(y)/dy-1)<1e-12);     % true if Y is linearly spaced

if ~evenx || ~eveny
  warning('(one of the) axis is/are not evenly spaced, but plots are made as if axis are linear')
end

% Select the channels in the data that match with the layout:
[seldat, sellay] = match_str(data.label, lay.label);
if isempty(seldat)
  error('labels in data and labels in layout do not match'); 
end

datavector = data.(cfg.zparam)(seldat,yidc,xidc);
chanX = lay.pos(sellay, 1);
chanY = lay.pos(sellay, 2);
chanWidth  = lay.width(sellay);
chanHeight = lay.height(sellay);
chanLabels = lay.label(sellay);
if ~isempty(cfg.maskparameter)
  maskvector = data.(cfg.maskparameter)(seldat,yidc,xidc);
end

% Get physical z-axis range (color axis):
if strcmp(cfg.zlim,'maxmin')
  zmin = min(datavector(:));
  zmax = max(datavector(:));
elseif strcmp(cfg.zlim,'absmax')
  zmin = -max(abs(datavector(:)));
  zmax = max(abs(datavector(:)));
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

hold on;

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

% Plot channels:
for k=1:length(seldat)
  % Get cdata:
  cdata = squeeze(datavector(k,:,:));
  if ~isempty(cfg.maskparameter)
    mdata = squeeze(maskvector(k,:,:));
  end
  
  % Get axes for this panel
  xas = (chanX(k) + linspace(0,1,size(cdata,2))*chanWidth(k)) - chanWidth(k)/2;
  yas = (chanY(k) + linspace(0,1,size(cdata,1))*chanHeight(k)) - chanHeight(k)/2;
  
  % Draw plot:
  h = imagesc(xas, yas, cdata, [zmin zmax]);
  
  % Mask Nan's and maskfield
  if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = double(mask);
    set(h,'AlphaData',mask, 'AlphaDataMapping', 'scaled');
    alim([0 1]);
  elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = mask .* mdata;
    mask = double(mask);
    set(h,'AlphaData',mask, 'AlphaDataMapping', 'scaled');
    alim([0 1]);
  elseif isequal(cfg.masknans,'no') && ~isempty(cfg.maskparameter)
    mask = mdata;
    mask = double(mask);
    set(h,'AlphaData',mask, 'AlphaDataMapping', 'scaled');
    alim([0 1]);
  end
  
% Draw box around plot
  if strcmp(cfg.box,'yes')
    xstep = xas(2) - xas(1); ystep = yas(2) - yas(1);
    xvalmin(1:length(yas)+2) = min(xas)-(0.5*xstep); xvalmax(1:length(yas)+2) = max(xas)+(0.5*xstep); yvalmin(1:length(xas)+2) = min(yas)-(0.5*ystep); yvalmax(1:length(xas)+2) = max(yas)+(0.5*ystep);
    xas2 = [xvalmin(1) xas xvalmax(1)]; yas2 = [yvalmin(1) yas yvalmax(1)];
    hold on
    plot([xas2 xvalmax xas2],[yvalmin yas2 yvalmax],'k');
    plot(xvalmin, yas2,'k');
  end

  % Draw channel labels:
  if strcmp(cfg.showlabels,'yes')
    text(chanX(k)-chanWidth(k)/2, chanY(k)+chanHeight(k)/2, sprintf(' %0s\n ', chanLabels{k}), 'Fontsize', cfg.fontsize);
  end
end

% write comment:
k = cellstrmatch('COMNT',lay.label);
if ~isempty(k)
  comment = cfg.comment;
  comment = sprintf('%0s\nxlim=[%.3g %.3g]', comment, xmin, xmax);
  comment = sprintf('%0s\nylim=[%.3g %.3g]', comment, ymin, ymax);
  comment = sprintf('%0s\nzlim=[%.3g %.3g]', comment, zmin, zmax);
  text(lay.pos(k,1), lay.pos(k,2), sprintf(comment), 'Fontsize', cfg.fontsize);
end

% plot scale:
k = cellstrmatch('SCALE',lay.label);
if ~isempty(k)
  % Get average cdata across channels:
  cdata = squeeze(mean(datavector, 1));
  
  % Get axes for this panel:
  xas = (lay.pos(k,1) + linspace(0,1,size(cdata,2))*lay.width(k));
  yas = (lay.pos(k,2) + linspace(0,1,size(cdata,1))*lay.height(k));
  
  % Draw plot:
  imagesc(xas, yas, cdata, [zmin zmax]);
end

% set colormap
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('singleplotTFR(): Colormap must be a n x 3 matrix'); end
  colormap(cfg.colormap);
end;

% plot colorbar:
if isfield(cfg, 'colorbar') && (strcmp(cfg.colorbar, 'yes'))
  colorbar;
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
    % add the channel information to the figure
    info       = guidata(gcf);
    info.x     = lay.pos(:,1);
    info.y     = lay.pos(:,2);
    info.label = lay.label;
    guidata(gcf, info);

    set(gcf, 'WindowButtonUpFcn',     {@select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
end

axis tight
axis off
if strcmp(cfg.box, 'yes')
  abc = axis;
  axis(abc + [-1 +1 -1 +1]*mean(abs(abc))/10)
end
orient landscape
hold off

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

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
% SUBFUNCTION which is called by select_channel in case cfg.cohrefchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_multiplotTFR(label, cfg, varargin)
cfg.cohrefchannel = label;
fprintf('selected cfg.cohrefchannel = ''%s''\n', cfg.cohrefchannel);
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
multiplotTFR(cfg, varargin{:});

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

