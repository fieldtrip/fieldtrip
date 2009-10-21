function [cfg] = singleplotTFR(cfg, data)

% singleplotTFR plots the time-frequency representations of power of a
% single channel or the average over multiple channels.
%
% Use as:
%   singleplotTFR(cfg,data)
%
% The data can be a time-frequency representation of power that was
% computed using the FREQANALYSIS function.
%
% The configuration can have the following parameters:
% cfg.xparam        = field to be plotted on x-axis, e.g. 'time' (default depends on data.dimord)
% cfg.yparam        = field to be plotted on y-axis, e.g. 'freq' (default depends on data.dimord)
% cfg.zparam        = field to be plotted on y-axis, e.g. 'powspcrtrm' (default depends on data.dimord)
% cfg.maskparameter = field in the data to be used for opacity masking of data
%                     (not possible for mean over multiple channels)
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.zlim          = 'maxmin','absmax' or [zmin zmax] (default = 'maxmin')
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see FREQBASELINE
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                     see CHANNELSELECTION for details
% cfg.cohrefchannel = name of reference channel for visualising coherence, can be 'gui'
% cfg.fontsize      = font size of title (default = 8)
% cfg.colormap      = any sized colormap, see COLORMAP
% cfg.colorbar      = 'yes', 'no' (default = 'yes')
% cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'no')
%                     In a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. Multiple areas
%                     can be selected by holding down the SHIFT key.
% cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
% cfg.masknans      = 'yes' or 'no' (default = 'yes')
%
% See also:
%   singleplotER, multiplotER, multiplotTFR, topoplotER, topoplotTFR.

% This function depends on FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype, documented

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% $Log: singleplotTFR.m,v $
% Revision 1.37  2009/07/14 13:52:16  roboos
% changed the interactive plotting: instead of using plotSelection it now uses the selection function from the new plotting module (select_range and select_channel) and uses a local subfunction to update the cfg and call the next figure
%
% Revision 1.36  2009/07/14 13:27:25  roboos
% consistent handling of cfg.renderer, default is to let matlab decide
%
% Revision 1.35  2009/06/17 13:44:52  roboos
% cleaned up help
%
% Revision 1.34  2009/05/12 18:47:23  roboos
% added general suppoprt for cfg.cohrefchannel and
% added handling of cfg.cohrefchannel='gui' for manual/interactive selection
%
% Revision 1.33  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.32  2008/12/16 14:59:12  sashae
% plot functions can now give cfg as output
% added checkconfig to start and end of function, configtracking possible
%
% Revision 1.31  2008/12/16 13:17:06  sashae
% replaced backward compatibility code by call to checkconfig
% removed obsolete subfunction pwrspctm2cohspctrm
%
% Revision 1.30  2008/11/28 16:33:58  sashae
% allow averaging over rpt/subj also for other fields than zparam=powspctrm (thanks to Jurrian)
%
% Revision 1.29  2008/10/28 14:24:05  ingnie
% fixed bug: data.time/freq should be data.(cfg.xparam/yparam)
%
% Revision 1.28  2008/10/27 12:00:02  ingnie
% change imagesc back to uimagesc, plotting with non linearly spaced axis possible now.
%
% Revision 1.27  2008/10/27 10:56:39  ingnie
% temporarily changed back uimagesc to imagesc, because adding of uimage files did not work
%
% Revision 1.26  2008/10/27 10:01:50  ingnie
% switched from using matlabs IMAGESC to UIMAGESC. Now also non linear spaced axis possible.
% Also simplified masking code, which made local variable -masking- unnecessary.
%
% Revision 1.25  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.24  2008/01/29 19:43:33  sashae
% added option for trial selection; plot functions now also accept data with
% repetitions (either trials or subjects), the avg is computed and plotted
% removed some old code
%
% Revision 1.23  2007/06/19 15:57:51  ingnie
% fixed bug in cfg.maskparameter, thanks to Saskia
%
% Revision 1.22  2007/06/19 14:01:54  ingnie
% added cfg.maskparameter, changed some white spaces
%
% Revision 1.21  2007/06/14 12:23:48  ingnie
% added cfg.colormap option
%
% Revision 1.20  2007/04/26 09:58:46  ingnie
% default masknans to 'yes'
%
% Revision 1.19  2007/04/25 17:25:06  ingnie
% added cfg.masknans option
%
% Revision 1.18  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.17  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.16  2007/01/09 10:41:43  roboos
% Added the option cfg.renderer, default is opengl. Interactive plotting
% on linux/VNC sometimes does not work, using cfg.renderer='painters'
% seems to fix it.
%
% Revision 1.15  2006/10/24 12:09:54  ingnie
% added colorbar yes/no option
%
% Revision 1.14  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.13  2006/06/07 10:24:34  ingnie
% added error when no channels are selected
%
% Revision 1.12  2006/06/06 16:24:06  ingnie
% replaced cfg.channelindex and cfg.channelname with cfg.channel for consistency
% with other functions
%
% Revision 1.11  2006/05/30 14:18:02  ingnie
% fixed bug that appeared when plotting single channel, updated documentation
%
% Revision 1.10  2006/05/09 17:33:57  ingnie
% fixed bug that appeared when cfg.channelindex is more than one channel
%
% Revision 1.9  2006/05/03 08:12:51  ingnie
% updated documentation
%
% Revision 1.8  2006/04/20 09:57:53  roboos
% changed formatting of the code from DOS into UNIX, i.e. removed the <CR>
%
% Revision 1.7  2006/04/07 23:45:23  chrhes
% changed clf to cla at beginning of function to allow use in subplots
%
% Revision 1.6  2006/03/22 18:56:57  jansch
% removed hard-coded selection of the powerspectrum in the to be plotted data
%
% Revision 1.5  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

cla

% For backward compatibility with old data structures:
data = checkdata(data);

% Set the defaults:
if ~isfield(cfg,'baseline'),        cfg.baseline = 'no';               end
if ~isfield(cfg,'baselinetype'),    cfg.baselinetype = 'absolute';     end
if ~isfield(cfg,'trials'),          cfg.trials = 'all';                end
if ~isfield(cfg,'xlim'),            cfg.xlim = 'maxmin';               end
if ~isfield(cfg,'ylim'),            cfg.ylim = 'maxmin';               end
if ~isfield(cfg,'zlim'),            cfg.zlim = 'maxmin';               end
if ~isfield(cfg,'fontsize'),        cfg.fontsize = 8;                  end
if ~isfield(cfg,'colorbar'),        cfg.colorbar = 'yes';              end
if ~isfield(cfg,'interactive'),     cfg.interactive = 'no';            end
if ~isfield(cfg,'renderer'),        cfg.renderer = [];                 end
if ~isfield(cfg,'masknans'),        cfg.masknans = 'yes';              end
if ~isfield(cfg,'maskparameter'),   cfg.maskparameter = [];            end

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

% Pick the channel(s)
if ~isfield(cfg,'channel')
  % set the default
  cfg.channel = 'all';
  % for backward compatibility
  cfg = checkconfig(cfg, 'renamed', {'channelindex',  'channel'});
  cfg = checkconfig(cfg, 'renamed', {'channelname',   'channel'});
end

% Check for unconverted coherence spectrum data
if (strcmp(cfg.zparam,'cohspctrm')) && (isfield(data, 'labelcmb'))
  % A reference channel is required:
  if ~isfield(cfg,'cohrefchannel'),
    error('no reference channel specified');
  end

  if strcmp(cfg.cohrefchannel, 'gui')
    % Open a single figure with the channel layout, the user can click on a reference channel
    h = clf;
    lay = prepare_layout(cfg, data);
    cfg.layout = lay;
    plot_lay(lay, 'box', false);
    title('Select the reference channel by clicking on it...');
    % add the channel information to the figure
    info       = guidata(h);
    info.x     = lay.pos(:,1);
    info.y     = lay.pos(:,2);
    info.label = lay.label;
    guidata(h, info);
    set(gcf, 'WindowButtonUpFcn', {@select_channel, 'callback', {@select_singleplotTFR, cfg, data}});
    return
  end

  % Convert 2-dimensional channel matrix to a single dimension:
  sel1                  = strmatch(cfg.cohrefchannel, data.labelcmb(:,2));
  sel2                  = strmatch(cfg.cohrefchannel, data.labelcmb(:,1));
  fprintf('selected %d channels for coherence\n', length(sel1)+length(sel2));
  data.cohspctrm = data.cohspctrm([sel1;sel2],:,:);
  data.label     = [data.labelcmb(sel1,1);data.labelcmb(sel2,2)];
  data.labelcmb  = data.labelcmb([sel1;sel2],:);
  data           = rmfield(data, 'labelcmb');
end

cfg.channel = channelselection(cfg.channel, data.label);
if isempty(cfg.channel)
  error('no channels selected');
else
  chansel = match_str(data.label, cfg.channel);
end

% cfg.maskparameter only possible for single channel
if length(chansel) > 1 && ~isempty(cfg.maskparameter)
  warning('no masking possible for average over multiple channels -> cfg.maskparameter cleared')
  cfg.maskparameter = [];
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

% Get TFR data averaged across selected channels, within the selected x/y-range:
dat = getsubfield(data, cfg.zparam);
TFR = squeeze(mean(dat(chansel,yidc,xidc), 1));
if ~isempty(cfg.maskparameter)
  mas = getsubfield(data, cfg.maskparameter);
  mdata = squeeze(mas(chansel,yidc,xidc));
end

% Get physical z-axis range (color axis):
if strcmp(cfg.zlim,'maxmin')
  zmin = min(TFR(:));
  zmax = max(TFR(:));
elseif strcmp(cfg.zlim,'absmax')
  zmin = -max(abs(TFR(:)));
  zmax = max(abs(TFR(:)));
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% test if X and Y are linearly spaced (to within 10^-12): % FROM UIMAGE
x = data.(cfg.xparam)(xidc);
y = data.(cfg.yparam)(yidc);
dx = min(diff(x));  % smallest interval for X
dy = min(diff(y));  % smallest interval for Y
evenx = all(abs(diff(x)/dx-1)<1e-12);     % true if X is linearly spaced
eveny = all(abs(diff(y)/dy-1)<1e-12);     % true if Y is linearly spaced

% masking only possible for evenly spaced axis
if strcmp(cfg.masknans, 'yes') && (~evenx || ~eveny)
  warning('(one of the) axis are not evenly spaced -> nans cannot be masked out ->  cfg.masknans is set to ''no'';')
  cfg.masknans = 'no';
end
if ~isempty(cfg.maskparameter) && (~evenx || ~eveny)
  warning('(one of the) axis are not evenly spaced -> no masking possible -> cfg.maskparameter cleared')
  cfg.maskparameter = [];
end

% Draw plot:
hold on;
h = uimagesc(data.(cfg.xparam)(xidc), data.(cfg.yparam)(yidc), TFR, [zmin,zmax]);
% Mask Nan's and maskfield
if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter)
  mask = ~isnan(TFR);
  mask = double(mask);
  set(h,'AlphaData',mask, 'AlphaDataMapping', 'scaled');
  alim([0 1]);
elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter)
  mask = ~isnan(TFR);
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
axis xy;

% set colormap
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('singleplotTFR(): Colormap must be a n x 3 matrix'); end
  colormap(cfg.colormap);
end;

if isequal(cfg.colorbar,'yes')
  colorbar;
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  set(gcf, 'WindowButtonUpFcn',     {@select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
end

% Create title text containing channel name(s) and channel number(s):
if length(chansel) == 1
  t = [char(cfg.channel) ' / ' num2str(chansel) ];
else
  t = sprintf('mean(%0s)', join(',', cfg.channel));
end
h = title(t,'fontsize', cfg.fontsize);

axis tight;
hold off;

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
function select_singleplotTFR(label, cfg, varargin)
cfg.cohrefchannel = label;
fprintf('selected cfg.cohrefchannel = ''%s''\n', cfg.cohrefchannel);
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
singleplotTFR(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotTFR(range, cfg, varargin)
cfg.comment = 'auto';
cfg.xlim = range(1:2);
cfg.ylim = range(3:4);
fprintf('selected cfg.xlim = [%f %f]\n', cfg.xlim(1), cfg.xlim(2));
fprintf('selected cfg.ylim = [%f %f]\n', cfg.ylim(1), cfg.ylim(2));
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
topoplotER(cfg, varargin{:});
