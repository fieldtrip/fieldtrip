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
% cfg.zparam        = field to be plotted on y-axis, e.g. 'powspcrtrm' (default depends on data.dimord)
% cfg.maskparameter = field in the data to be used for masking of data
%                     (not possible for mean over multiple channels, or when input contains multiple subjects
%                     or trials)
% cfg.maskstyle     = style used to mask nans, 'opacity' or 'saturation' (default = 'opacity')
%                     use 'saturation' when saving to vector-format (like *.eps) to avoid all sorts of image-problems
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

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

cla

% For backward compatibility with old data structures:
data = ft_checkdata(data);

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'renamedval',  {'zlim',  'absmax',  'maxabs'});

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
if ~isfield(cfg,'maskstyle'),       cfg.maskstyle = 'opacity';         end
if ~isfield(cfg, 'matrixside'),     cfg.matrixside = 'feedforward';    end

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

% Pick the channel(s)
if ~isfield(cfg,'channel')
  % set the default
  cfg.channel = 'all';
  % for backward compatibility
  cfg = checkconfig(cfg, 'renamed', {'channelindex',  'channel'});
  cfg = checkconfig(cfg, 'renamed', {'channelname',   'channel'});
end

cfg.channel = ft_channelselection(cfg.channel, data.label);
if isempty(cfg.channel)
  error('no channels selected');
else
  chansel = match_str(data.label, cfg.channel);
end

% Check for unconverted coherence spectrum data or any other bivariate metric:
dimtok  = tokenize(data.dimord, '_');
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;
if (strcmp(cfg.zparam,'cohspctrm')) && (isfield(data, 'labelcmb')) || ...
    (isfull && isfield(data, cfg.zparam)),

  % A reference channel is required:
  if ~isfield(cfg,'cohrefchannel'),
    error('no reference channel specified');
  end

  if strcmp(cfg.cohrefchannel, 'gui')
    % Open a single figure with the channel layout, the user can click on a reference channel
    h = clf;
    lay = ft_prepare_layout(cfg, data);
    cfg.layout = lay;
    ft_plot_lay(lay, 'box', false);
    title('Select the reference channel by clicking on it...');
    % add the channel information to the figure
    info       = guidata(h);
    info.x     = lay.pos(:,1);
    info.y     = lay.pos(:,2);
    info.label = lay.label;
    guidata(h, info);
    set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'callback', {@select_singleplotTFR, cfg, data}});
    return
  end


  if ~isfull,
    % only works explicitly with coherence FIXME
    % Convert 2-dimensional channel matrix to a single dimension:
    sel1           = strmatch(cfg.cohrefchannel, data.labelcmb(:,2));
    sel2           = strmatch(cfg.cohrefchannel, data.labelcmb(:,1));
    fprintf('selected %d channels for coherence\n', length(sel1)+length(sel2));
    data.cohspctrm = data.cohspctrm([sel1;sel2],:,:);
    data.label     = [data.labelcmb(sel1,1);data.labelcmb(sel2,2)];
    data.labelcmb  = data.labelcmb([sel1;sel2],:);
    data           = rmfield(data, 'labelcmb');
  else
    % general solution
    
    sel               = match_str(data.label, cfg.cohrefchannel);
    siz               = [size(data.(cfg.zparam)) 1];
    if strcmp(cfg.matrixside, 'feedback')
      %data.(cfg.zparam) = reshape(mean(data.(cfg.zparam)(:,sel,:),2),[siz(1) 1 siz(3:end)]);
      %sel1 = 1:siz(1);
      sel1 = chansel;
      sel2 = sel;
      meandir = 2;
    elseif strcmp(cfg.matrixside, 'feedforward')
      %data.(cfg.zparam) = reshape(mean(data.(cfg.zparam)(sel,:,:),1),[siz(1) 1 siz(3:end)]);
      sel1 = sel;
      %sel2 = 1:siz(1);
      sel2 = chansel;
      meandir = 1;
    elseif strcmp(cfg.matrixside, 'ff-fd')
      %FIXME don't know how to handle this
      data.(cfg.zparam) = reshape(mean(data.(cfg.zparam)(sel,:,:),1),[siz(1) 1 siz(3:end)]) - reshape(mean(data.(cfg.zparam)(:,sel,:),2),[siz(1) 1 siz(3:end)]);
    elseif strcmp(cfg.matrixside, 'fd-ff')
      data.(cfg.zparam) = reshape(mean(data.(cfg.zparam)(:,sel,:),2),[siz(1) 1 siz(3:end)]) - reshape(mean(data.(cfg.zparam)(sel,:,:),1),[siz(1) 1 siz(3:end)]);
    end
  end
end

% cfg.maskparameter only possible for single channel
if length(chansel) > 1 && ~isempty(cfg.maskparameter)
  warning('no masking possible for average over multiple channels -> cfg.maskparameter cleared')
  cfg.maskparameter = [];
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

% Get TFR data averaged across selected channels, within the selected x/y-range:
if ~isfull
  dat = getsubfield(data, cfg.zparam);
  TFR = squeeze(mean(dat(chansel,yidc,xidc), 1));
  if ~isempty(cfg.maskparameter)
    mas = getsubfield(data, cfg.maskparameter);
    mdata = squeeze(mas(chansel,yidc,xidc));
  end
elseif isfull,
  %siz = size(data.(cfg.zparam));
  %siz(meandir) = [];
  TFR = squeeze(mean(mean(data.(cfg.zparam)(sel1,sel2,yidc,xidc),1),2));
end

% Get physical z-axis range (color axis):
if strcmp(cfg.zlim,'maxmin')
  zmin = min(TFR(:));
  zmax = max(TFR(:));
elseif strcmp(cfg.zlim,'maxabs')
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

% set colormap
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('singleplotTFR(): Colormap must be a n x 3 matrix'); end
  set(gcf,'colormap',cfg.colormap);
end;

% Draw plot (and mask NaN's if requested):
if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter)
  mask = ~isnan(TFR);
  mask = double(mask);
  ft_plot_matrix(data.(cfg.xparam)(xidc),data.(cfg.yparam)(yidc), TFR, 'clim',[zmin,zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter)
  mask = ~isnan(TFR);
  mask = mask .* mdata;
  mask = double(mask);
  ft_plot_matrix(data.(cfg.xparam)(xidc),data.(cfg.yparam)(yidc), TFR, 'clim',[zmin,zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
elseif isequal(cfg.masknans,'no') && ~isempty(cfg.maskparameter)
  mask = mdata;
  mask = double(mask);
  ft_plot_matrix(data.(cfg.xparam)(xidc),data.(cfg.yparam)(yidc), TFR, 'clim',[zmin,zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
else
  ft_plot_matrix(data.(cfg.xparam)(xidc),data.(cfg.yparam)(yidc), TFR, 'clim',[zmin,zmax],'tag','cip')
end
hold on
axis xy;



if isequal(cfg.colorbar,'yes')
  colorbar;
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
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

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end


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
