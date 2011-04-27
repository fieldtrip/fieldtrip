function [cfg] = ft_singleplotTFR(cfg, data)

% ft_singleplotTFR plots the time-frequency representations of power of a
% single channel or the average over multiple channels.
%
% Use as:
%   ft_singleplotTFR(cfg,data)
%
% The data can be a time-frequency representation of power that was
% computed using the FT_FREQANALYSIS function.
%
% The configuration can have the following parameters:
% cfg.xparam        = field to be plotted on x-axis, e.g. 'time' (default depends on data.dimord)
% cfg.yparam        = field to be plotted on y-axis, e.g. 'freq' (default depends on data.dimord)
% cfg.zparam        = field to be plotted on z-axis, e.g. 'powspcrtrm' (default depends on data.dimord)
% cfg.maskparameter = field in the data to be used for masking of data
%                     (not possible for mean over multiple channels, or when input contains multiple subjects
%                     or trials)
% cfg.maskstyle     = style used to mask nans, 'opacity' or 'saturation' (default = 'opacity')
%                     use 'saturation' when saving to vector-format (like *.eps) to avoid all sorts of image-problems
% cfg.maskalpha     = alpha value used for masking areas dictated by cfg.maskparameter (0 - 1, default = 1)
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.zlim          = 'maxmin','maxabs' or [zmin zmax] (default = 'maxmin')
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see FT_FREQBASELINE
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
% cfg.cohrefchannel = name of reference channel for visualising coherence, can be 'gui'
% cfg.fontsize      = font size of title (default = 8)
% cfg.hotkeys          = enables hotkeys (up/down arrows) for dynamic colorbar adjustment
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
%   FT_SINGLEPLOTER, FT_MULTIPLOTER, FT_MULTIPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR

% This function depends on FT_FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype, documented

% Copyright (C) 2005-2006, F.C. Donders Centre
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
cfg = ft_checkconfig(cfg, 'unused',      {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'zlim',  'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'renamed',     {'channelindex',  'channel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'channelname',   'channel'});

cla

% Set the defaults:
cfg.baseline      = ft_getopt(cfg, 'baseline',     'no');
cfg.baselinetype  = ft_getopt(cfg, 'baselinetype', 'absolute');
cfg.trials        = ft_getopt(cfg, 'trials',       'all');
cfg.xlim          = ft_getopt(cfg, 'xlim',         'maxmin'); 
cfg.ylim          = ft_getopt(cfg, 'ylim',         'maxmin');
cfg.zlim          = ft_getopt(cfg, 'zlim',         'maxmin');
cfg.fontsize      = ft_getopt(cfg, 'fontsize',     8);
cfg.colorbar      = ft_getopt(cfg, 'colorbar',     'yes');
cfg.interactive   = ft_getopt(cfg, 'interactive',  'no');
cfg.hotkeys       = ft_getopt(cfg, 'hotkeys',      'no');
cfg.renderer      = ft_getopt(cfg, 'renderer',     []);
cfg.maskalpha     = ft_getopt(cfg, 'maskalpha',     1);
cfg.maskparameter = ft_getopt(cfg, 'maskparameter',[]);
cfg.maskstyle     = ft_getopt(cfg, 'maskstyle',    'opacity');
cfg.channel       = ft_getopt(cfg, 'channel',      'all');
cfg.masknans      = ft_getopt(cfg, 'masknans',     'yes');
cfg.matrixside    = ft_getopt(cfg, 'matrixside',   '');

% for backward compatibility with old data structures
data   = ft_checkdata(data, 'datatype', 'freq');
dimord = data.dimord;
dimtok = tokenize(dimord, '_');

% Set x/y/zparam defaults
if ~sum(ismember(dimtok, 'time'))
  error('input data needs a time dimension');
else
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
end

if isfield(cfg, 'channel') && isfield(data, 'label')
  cfg.channel = ft_channelselection(cfg.channel, data.label);
elseif isfield(cfg, 'channel') && isfield(data, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
end

if isempty(cfg.channel)
  error('no channels selected');
end

% check whether rpt/subj is present and remove if necessary and whether
hasrpt = sum(ismember(dimtok, {'rpt' 'subj'}));
if hasrpt,
  % this also deals with fourier-spectra in the input
  % or with multiple subjects in a frequency domain stat-structure
  % on the fly computation of coherence spectrum is not supported
  if isfield(data, 'crsspctrm'),
    data = rmfield(data, 'crsspctrm'); 
  end
  
  tmpcfg           = [];
  tmpcfg.trials    = cfg.trials;
  tmpcfg.jackknife = 'no';
  if isfield(cfg, 'zparam') && ~strcmp(cfg.zparam,'powspctrm')
    % freqdesctiptives will only work on the powspctrm field
    % hence a temporary copy of the data is needed
    tempdata.dimord    = data.dimord;
    tempdata.freq      = data.freq;
    tempdata.label     = data.label;
    tempdata.powspctrm = data.(cfg.zparam);
    tempdata.cfg       = data.cfg;
    tempdata           = ft_freqdescriptives(tmpcfg, tempdata);
    data.(cfg.zparam)  = tempdata.powspctrm;
    clear tempdata
  else
    data = ft_freqdescriptives(tmpcfg, data);
  end
  dimord = data.dimord;
  dimtok = tokenize(dimord, '_');
end % if hasrpt

% Handle the bivariate case

% Check for bivariate metric with 'chan_chan' in the dimord
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;

% Check for bivariate metric with a labelcmb
haslabelcmb = isfield(data, 'labelcmb');

if (isfull || haslabelcmb) && isfield(data, cfg.zparam)
  % A reference channel is required:
  if ~isfield(cfg, 'cohrefchannel')
    error('no reference channel is specified');
  end
  
  % check for cohrefchannel being part of selection
  if ~strcmp(cfg.cohrefchannel,'gui')
    if (isfull      && ~any(ismember(data.label, cfg.cohrefchannel))) || ...
       (haslabelcmb && ~any(ismember(data.labelcmb(:), cfg.cohrefchannel)))
      error('cfg.cohrefchannel is a not present in the (selected) channels)')
    end
  end
  
  % Interactively select the reference channel
  if strcmp(cfg.cohrefchannel, 'gui')
    error('coh.refchannel = ''gui'' is not supported at the moment for ft_singleplotTFR');
%     
%     % Open a single figure with the channel layout, the user can click on a reference channel
%     h = clf;
%     ft_plot_lay(lay, 'box', false);
%     title('Select the reference channel by dragging a selection window, more than 1 channel can be selected...');
%     % add the channel information to the figure
%     info       = guidata(gcf);
%     info.x     = lay.pos(:,1);
%     info.y     = lay.pos(:,2);
%     info.label = lay.label;
%     guidata(h, info);
%     %set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'callback', {@select_topoplotER, cfg, data}});
%     set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
%     set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
%     set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
%     return
  end
  
  if ~isfull,
    % Convert 2-dimensional channel matrix to a single dimension:
    if isempty(cfg.matrixside)
      sel1 = strmatch(cfg.cohrefchannel, data.labelcmb(:,2), 'exact');
      sel2 = strmatch(cfg.cohrefchannel, data.labelcmb(:,1), 'exact');
    elseif strcmp(cfg.matrixside, 'feedforward')
      sel1 = [];
      sel2 = strmatch(cfg.cohrefchannel, data.labelcmb(:,1), 'exact');
    elseif strcmp(cfg.matrixside, 'feedback')
      sel1 = strmatch(cfg.cohrefchannel, data.labelcmb(:,2), 'exact');
      sel2 = [];
    end
    fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.zparam);
    data.(cfg.zparam) = data.(cfg.zparam)([sel1;sel2],:,:);
    data.label     = [data.labelcmb(sel1,1);data.labelcmb(sel2,2)];
    data.labelcmb  = data.labelcmb([sel1;sel2],:);
    data           = rmfield(data, 'labelcmb');
  else
    % General case
    sel               = match_str(data.label, cfg.cohrefchannel);
    siz               = [size(data.(cfg.zparam)) 1];
    if strcmp(cfg.matrixside, 'feedback') || isempty(cfg.matrixside)
      %FIXME the interpretation of 'feedback' and 'feedforward' depend on
      %the definition in the bivariate representation of the data
      %data.(cfg.zparam) = reshape(mean(data.(cfg.zparam)(:,sel,:),2),[siz(1) 1 siz(3:end)]);
      sel1 = 1:siz(1);
      sel2 = sel;
      meandir = 2;
    elseif strcmp(cfg.matrixside, 'feedforward')
      %data.(cfg.zparam) = reshape(mean(data.(cfg.zparam)(sel,:,:),1),[siz(1) 1 siz(3:end)]);
      sel1 = sel;
      sel2 = 1:siz(1);
      meandir = 1;

    elseif strcmp(cfg.matrixside, 'ff-fd')
      error('cfg.matrixside = ''ff-fd'' is not supported anymore, you have to manually subtract the two before the call to ft_topoplotER');
    elseif strcmp(cfg.matrixside, 'fd-ff')
      error('cfg.matrixside = ''fd-ff'' is not supported anymore, you have to manually subtract the two before the call to ft_topoplotER');
    end %if matrixside
  end %if ~isfull
end %handle the bivariate data

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

% Replace value with the index of the nearest bin
if ~isempty(cfg.xparam)
  xmin = nearest(data.(cfg.xparam), xmin);
  xmax = nearest(data.(cfg.xparam), xmax);
end

% Get physical y-axis range:
if strcmp(cfg.ylim,'maxmin')
  ymin = min(data.(cfg.yparam));
  ymax = max(data.(cfg.yparam));
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Replace value with the index of the nearest bin
if ~isempty(cfg.yparam)
  ymin = nearest(data.(cfg.yparam), ymin);
  ymax = nearest(data.(cfg.yparam), ymax);
end

% test if X and Y are linearly spaced (to within 10^-12): % FROM UIMAGE
x = data.(cfg.xparam)(xmin:xmax);
y = data.(cfg.yparam)(ymin:ymax);
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

% perform channel selection
selchannel = ft_channelselection(cfg.channel, data.label);
sellab     = match_str(data.label, selchannel);

% cfg.maskparameter only possible for single channel
if length(sellab) > 1 && ~isempty(cfg.maskparameter)
  warning('no masking possible for average over multiple channels -> cfg.maskparameter cleared')
  cfg.maskparameter = [];
end

dat = data.(cfg.zparam);
if isfull
  dat = dat(sel1, sel2, ymin:ymax, xmin:xmax);
  dat = nanmean(dat, meandir);
  siz = size(dat);
  dat = reshape(dat, [max(siz(1:2)) siz(3) siz(4)]);
  dat = dat(sellab, :, :);
elseif haslabelcmb
  dat = dat(sellab, ymin:ymax, xmin:xmax);
else
  dat = dat(sellab, ymin:ymax, xmin:xmax);
end

if ~isempty(cfg.maskparameter)
  mask = data.(cfg.maskparameter);
  if isfull && cfg.maskalpha == 1
    mask = mask(sel1, sel2, ymin:ymax, xmin:xmax);
    mask = nanmean(mask, meandir);
    siz  = size(mask);
    mask = reshape(mask, [max(siz(1:2)) siz(3) siz(4)]);
    mask = reshape(mask(sellab, :, :), [siz(3) siz(4)]);
  elseif haslabelcmb && cfg.maskalpha == 1
    mask = mask(sellab, ymin:ymax, xmin:xmax);
  elseif cfg.maskalpha == 1
    mask = mask(sellab, ymin:ymax, xmin:xmax);
  elseif isfull && cfg.maskalpha ~= 1 %% check me
    maskl = mask(sel1, sel2, ymin:ymax, xmin:xmax);
    maskl = nanmean(maskl, meandir);
    siz  = size(maskl);
    maskl = reshape(maskl, [max(siz(1:2)) siz(3) siz(4)]);
    maskl = squeeze(reshape(maskl(sellab, :, :), [siz(3) siz(4)]));
    mask = zeros(size(maskl));
    mask(maskl) = 1;
    mask(~maskl) = cfg.maskalpha;
  elseif haslabelcmb && cfg.maskalpha ~= 1
    maskl = squeeze(mask(sellab, ymin:ymax, xmin:xmax));
    mask = zeros(size(maskl));
    mask(maskl) = 1;
    mask(~maskl) = cfg.maskalpha;
  elseif cfg.maskalpha ~= 1
    maskl = squeeze(mask(sellab, ymin:ymax, xmin:xmax));
    mask = zeros(size(maskl));
    mask(maskl) = 1;
    mask(~maskl) = cfg.maskalpha;
  end
end
siz        = size(dat);
datamatrix = reshape(mean(dat, 1), siz(2:end));
xvector    = data.(cfg.xparam)(xmin:xmax);
yvector    = data.(cfg.yparam)(ymin:ymax);

% Get physical z-axis range (color axis):
if strcmp(cfg.zlim,'maxmin')
  zmin = min(datamatrix(:));
  zmax = max(datamatrix(:));
elseif strcmp(cfg.zlim,'maxabs')
  zmin = -max(abs(datamatrix(:)));
  zmax = max(abs(datamatrix(:)));
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% set colormap
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('singleplotTFR(): Colormap must be a n x 3 matrix'); end
  set(gcf,'colormap',cfg.colormap);
end;

% Draw plot (and mask NaN's if requested):
if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter)
  nans_mask = ~isnan(datamatrix);
  mask = double(nans_mask);
  ft_plot_matrix(xvector, yvector, datamatrix, 'clim',[zmin,zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter)
  nans_mask = ~isnan(datamatrix);
  mask = mask .* nans_mask;
  mask = double(mask);
  ft_plot_matrix(xvector, yvector, datamatrix, 'clim',[zmin,zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
elseif isequal(cfg.masknans,'no') && ~isempty(cfg.maskparameter)
  mask = double(mask);
  ft_plot_matrix(xvector, yvector, datamatrix, 'clim',[zmin,zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
else
  ft_plot_matrix(xvector, yvector, datamatrix, 'clim',[zmin,zmax],'tag','cip')
end
hold on
axis xy;
set(gca,'Color','k')

if isequal(cfg.colorbar,'yes')
  colorbar;
end

% Set adjust color axis
if strcmp('yes',cfg.hotkeys)
  %  Attach data and cfg to figure and attach a key listener to the figure
  set(gcf, 'KeyPressFcn', {@key_sub, zmin, zmax})
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
end

% Create title text containing channel name(s) and channel number(s):
if length(sellab) == 1
  t = [char(cfg.channel) ' / ' num2str(sellab) ];
else
  t = sprintf('mean(%0s)', join(',', cfg.channel));
end
h = title(t,'fontsize', cfg.fontsize);

axis tight;
hold off;

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
function select_singleplotTFR(label, cfg, varargin)
cfg.cohrefchannel = label;
fprintf('selected cfg.cohrefchannel = ''%s''\n', cfg.cohrefchannel);
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_singleplotTFR(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotTFR(range, cfg, varargin)
cfg.comment = 'auto';
cfg.xlim = range(1:2);
cfg.ylim = range(3:4);
% compatibility fix for new ft_topoplotER/TFR cfg options
if isfield(cfg,'showlabels') && strcmp(cfg.showlabels,'yes')
  cfg = rmfield(cfg,'showlabels');
  cfg.marker = 'labels';
elseif isfield(cfg,'showlabels') && strcmp(cfg.showlabels,'no')
  cfg = rmfield(cfg,'showlabels');
  cfg.marker = 'on';
end
fprintf('selected cfg.xlim = [%f %f]\n', cfg.xlim(1), cfg.xlim(2));
fprintf('selected cfg.ylim = [%f %f]\n', cfg.ylim(1), cfg.ylim(2));
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_topoplotER(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata, varargin)
incr = (max(caxis)-min(caxis)) /10;
% symmetrically scale color bar down by 10 percent
if strcmp(eventdata.Key,'uparrow')
  caxis([min(caxis)-incr max(caxis)+incr]);
% symmetrically scale color bar up by 10 percent
elseif strcmp(eventdata.Key,'downarrow')
  caxis([min(caxis)+incr max(caxis)-incr]);
% resort to minmax of data for colorbar
elseif strcmp(eventdata.Key,'m')
  caxis([varargin{1} varargin{2}]);
end