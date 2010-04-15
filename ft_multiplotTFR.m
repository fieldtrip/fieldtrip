function [cfg] = ft_multiplotTFR(cfg, data)

% ft_multiplotTFR plots time-frequency representations of power or coherence in a 
% topographical layout. The plots of the indivual sensors are arranged according 
% to their location specified in the layout.
%
% Use as:
%   ft_multiplotTFR(cfg, data)
%
% The data can be a time-frequency representation of power or coherence that 
% was computed using the FT_FREQANALYSIS or FT_FREQDESCRIPTIVES functions.
%
% The configuration can have the following parameters:
% cfg.xparam           = field to be plotted on x-axis (default depends on data.dimord)
%                        'time'
% cfg.yparam           = field to be plotted on y-axis (default depends on data.dimord)
%                        'freq'
% cfg.zparam           = field to be represented as color (default depends on data.dimord)
%                        'powspctrm' or 'cohspctrm' 
% cfg.maskparameter    = field in the data to be used for opacity masking of data
% cfg.maskstyle        = style used to mask nans, 'opacity' or 'saturation' (default = 'opacity')
%                        use 'saturation' when saving to vector-format (like *.eps) to avoid all 
%                        sorts of image-problems (currently only possible with a white backgroud)
% cfg.xlim             = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim             = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.zlim             = 'maxmin','maxabs' or [zmin zmax] (default = 'maxmin')
% cfg.cohrefchannel    = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline         = 'yes','no' or [time1 time2] (default = 'no'), see FT_FREQBASELINE
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
%  - you can provide a pre-computed layout structure (see ft_prepare_layout)
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
%   ft_multiplotER, ft_singleplotER, ft_singleplotTFR, ft_topoplotER, ft_topoplotTFR,
%   ft_prepare_layout

% Undocumented local options:
% cfg.channel
% cfg.layoutname
% cfg.xparam
% cfg.zparam
%
% This function depends on FT_FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype, documented

% Copyright (C) 2003-2006, Ole Jensen
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

clf

% for backward compatibility with old data structures
data = checkdata(data);

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'renamedval',  {'zlim',  'absmax',  'maxabs'});

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
if ~isfield(cfg,'maskstyle'),       cfg.maskstyle = 'opacity';         end
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
    tempdata           = ft_freqdescriptives(tmpcfg, tempdata);
    data.(cfg.zparam)  = tempdata.powspctrm;
    clear tempdata
  else
    data = ft_freqdescriptives(tmpcfg, data);
  end
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
end

% Old style coherence plotting with cohtargetchannel is no longer supported:
cfg = checkconfig(cfg, 'unused',  {'cohtargetchannel'});

% Read or create the layout that will be used for plotting:
lay = ft_prepare_layout(cfg, data);
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
  data = ft_freqbaseline(cfg, data);
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
elseif strcmp(cfg.zlim,'maxabs')
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

% set colormap
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('singleplotTFR(): Colormap must be a n x 3 matrix'); end
  set(gcf,'colormap',cfg.colormap);
end;

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
  
  % Draw plot (and mask Nan's with maskfield if requested)
  if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = double(mask);
    plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = mask .* mdata;
    mask = double(mask);
    plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  elseif isequal(cfg.masknans,'no') && ~isempty(cfg.maskparameter)
    mask = mdata;
    mask = double(mask);
    plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  else
    plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip')
  end
  % Currently the handle isn't being used below, this is here for possible use in the future
  h = findobj('tag','cip');
  
  % Draw box around plot
  if strcmp(cfg.box,'yes')
    xstep = xas(2) - xas(1); ystep = yas(2) - yas(1);
    xvalmin(1:length(yas)+2) = min(xas)-(0.5*xstep); xvalmax(1:length(yas)+2) = max(xas)+(0.5*xstep); yvalmin(1:length(xas)+2) = min(yas)-(0.5*ystep); yvalmax(1:length(xas)+2) = max(yas)+(0.5*ystep);
    xas2 = [xvalmin(1) xas xvalmax(1)]; yas2 = [yvalmin(1) yas yvalmax(1)];
    plot_box([min(xas2) max(xas2) min(yas2) max(yas2)])
  end

  % Draw channel labels:
  if strcmp(cfg.showlabels,'yes')
    plot_text(chanX(k)-chanWidth(k)/2, chanY(k)+chanHeight(k)/2, sprintf(' %0s\n ', chanLabels{k}), 'Fontsize', cfg.fontsize);
  end
end

% write comment:
k = cellstrmatch('COMNT',lay.label);
if ~isempty(k)
  comment = cfg.comment;
  comment = sprintf('%0s\nxlim=[%.3g %.3g]', comment, xmin, xmax);
  comment = sprintf('%0s\nylim=[%.3g %.3g]', comment, ymin, ymax);
  comment = sprintf('%0s\nzlim=[%.3g %.3g]', comment, zmin, zmax);
  plot_text(lay.pos(k,1), lay.pos(k,2), sprintf(comment), 'Fontsize', cfg.fontsize);
end

% plot scale:
k = cellstrmatch('SCALE',lay.label);
if ~isempty(k)
  % Get average cdata across channels:
  cdata = squeeze(mean(datavector, 1));
  
  % Get axes for this panel:
  xas = (lay.pos(k,1) + linspace(0,1,size(cdata,2))*lay.width(k));
  yas = (lay.pos(k,2) + linspace(0,1,size(cdata,1))*lay.height(k));
  
  % Draw plot (and mask Nan's with maskfield if requested)
  if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = double(mask);
    plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = mask .* mdata;
    mask = double(mask);
    plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  elseif isequal(cfg.masknans,'no') && ~isempty(cfg.maskparameter)
    mask = mdata;
    mask = double(mask);
    plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  else
    plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip')
  end
  % Currently the handle isn't being used below, this is here for possible use in the future
  h = findobj('tag','cip');

end


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

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

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
ft_multiplotTFR(cfg, varargin{:});

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
  ft_singleplotTFR(cfg, varargin{:});
end

