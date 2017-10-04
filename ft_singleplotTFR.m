function [cfg] = ft_singleplotTFR(cfg, data)

% FT_SINGLEPLOTTFR plots the time-frequency representation of power of a
% single channel or the average over multiple channels.
%
% Use as
%   ft_singleplotTFR(cfg,data)
%
% The input freq structure should be a a time-frequency representation of
% power or coherence that was computed using the FT_FREQANALYSIS function.
%
% The configuration can have the following parameters:
%   cfg.parameter      = field to be plotted on z-axis, e.g. 'powspcrtrm' (default depends on data.dimord)
%   cfg.maskparameter  = field in the data to be used for masking of data
%                        (not possible for mean over multiple channels, or when input contains multiple subjects
%                        or trials)
%   cfg.maskstyle      = style used to masking, 'opacity', 'saturation', 'outline' or 'colormix' (default = 'opacity')
%                        use 'saturation' or 'outline' when saving to vector-format (like *.eps) to avoid all sorts of image-problems
%   cfg.maskalpha      = alpha value between 0 (transparant) and 1 (opaque) used for masking areas dictated by cfg.maskparameter (default = 1)
%   cfg.masknans       = 'yes' or 'no' (default = 'yes')
%   cfg.xlim           = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim           = 'maxmin' or [ymin ymax] (default = 'maxmin')
%   cfg.zlim           = plotting limits for color dimension, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.baseline       = 'yes', 'no' or [time1 time2] (default = 'no'), see FT_FREQBASELINE
%   cfg.baselinetype   = 'absolute', 'relative', 'relchange' or 'db' (default = 'absolute')
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                        see FT_CHANNELSELECTION for details
%   cfg.title          = string, title of plot
%   cfg.refchannel     = name of reference channel for visualising connectivity, can be 'gui'
%   cfg.fontsize       = font size of title (default = 8)
%   cfg.hotkeys        = enables hotkeys (leftarrow/rightarrow/uparrow/downarrow/pageup/pagedown/m) for dynamic zoom and translation (ctrl+) of the axes and color limits
%   cfg.colormap       = any sized colormap, see COLORMAP
%   cfg.colorbar       = 'yes', 'no' (default = 'yes')
%   cfg.interactive    = Interactive plot 'yes' or 'no' (default = 'yes')
%                        In a interactive plot you can select areas and produce a new
%                        interactive plot when a selected area is clicked. Multiple areas
%                        can be selected by holding down the SHIFT key.
%   cfg.renderer       = 'painters', 'zbuffer', ' opengl' or 'none' (default = [])
%   cfg.directionality = '', 'inflow' or 'outflow' specifies for
%                       connectivity measures whether the inflow into a
%                       node, or the outflow from a node is plotted. The
%                       (default) behavior of this option depends on the dimor
%                       of the input data (see below).
%
% For the plotting of directional connectivity data the cfg.directionality
% option determines what is plotted. The default value and the supported
% functionality depend on the dimord of the input data. If the input data
% is of dimord 'chan_chan_XXX', the value of directionality determines
% whether, given the reference channel(s), the columns (inflow), or rows
% (outflow) are selected for plotting. In this situation the default is
% 'inflow'. Note that for undirected measures, inflow and outflow should
% give the same output. If the input data is of dimord 'chancmb_XXX', the
% value of directionality determines whether the rows in data.labelcmb are
% selected. With 'inflow' the rows are selected if the refchannel(s) occur in
% the right column, with 'outflow' the rows are selected if the
% refchannel(s) occur in the left column of the labelcmb-field. Default in
% this case is '', which means that all rows are selected in which the
% refchannel(s) occur. This is to robustly support linearly indexed
% undirected connectivity metrics. In the situation where undirected
% connectivity measures are linearly indexed, specifying 'inflow' or
% 'outflow' can result in unexpected behavior.
%
% See also FT_SINGLEPLOTER, FT_MULTIPLOTER, FT_MULTIPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR

% Copyright (C) 2005-2017, F.C. Donders Centre
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPERS NOTE: This code is organized in a similar fashion for multiplot/singleplot/topoplot
% and for ER/TFR and should remain consistent over those 6 functions.
% Section 1: general cfg handling that is independent from the data
% Section 2: data handling, this also includes converting bivariate (chan_chan and chancmb) into univariate data
% Section 3: select the data to be plotted and determine min/max range
% Section 4: do the actual plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 1: general cfg handling that is independent from the data

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'freq');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'unused',      {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'matrixside',     'directionality'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'zlim', 'absmax', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'directionality', 'feedback',    'inflow'});
cfg = ft_checkconfig(cfg, 'renamed',     {'channelindex',   'channel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'channelname',    'channel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'cohrefchannel',  'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed',	   {'zparam',         'parameter'});
cfg = ft_checkconfig(cfg, 'deprecated',  {'xparam',         'yparam'});

% Set the defaults
cfg.baseline       = ft_getopt(cfg, 'baseline',      'no');
cfg.baselinetype   = ft_getopt(cfg, 'baselinetype',  'absolute');
cfg.trials         = ft_getopt(cfg, 'trials',        'all', 1);
cfg.xlim           = ft_getopt(cfg, 'xlim',          'maxmin');
cfg.ylim           = ft_getopt(cfg, 'ylim',          'maxmin');
cfg.zlim           = ft_getopt(cfg, 'zlim',          'maxmin');
cfg.fontsize       = ft_getopt(cfg, 'fontsize',       8);
cfg.colorbar       = ft_getopt(cfg, 'colorbar',      'yes');
cfg.interactive    = ft_getopt(cfg, 'interactive',   'yes');
cfg.hotkeys        = ft_getopt(cfg, 'hotkeys',       'yes');
cfg.renderer       = ft_getopt(cfg, 'renderer',       []);
cfg.maskalpha      = ft_getopt(cfg, 'maskalpha',      1);
cfg.maskparameter  = ft_getopt(cfg, 'maskparameter',  []);
cfg.maskstyle      = ft_getopt(cfg, 'maskstyle',     'opacity');
cfg.channel        = ft_getopt(cfg, 'channel',       'all');
cfg.title          = ft_getopt(cfg, 'title',          []);
cfg.masknans       = ft_getopt(cfg, 'masknans',      'yes');
cfg.directionality = ft_getopt(cfg, 'directionality', []);
cfg.figurename     = ft_getopt(cfg, 'figurename',     []);
cfg.parameter      = ft_getopt(cfg, 'parameter',     'powspctrm');

% this is needed for the figure title and correct labeling of graphcolor later on
if nargin>1
  if isfield(cfg, 'dataname')
    if iscell(cfg.dataname)
      dataname = cfg.dataname{1};
    else
      dataname = cfg.dataname;
    end
  else
    if ~isempty(inputname(2))
      dataname = inputname(2);
    else
      dataname = ['data' num2str(1, '%02d')];
    end
  end
else  % data provided through cfg.inputfile
  dataname = cfg.inputfile;
end


%% Section 2: data handling, this also includes converting bivariate (chan_chan and chancmb) into univariate data

hastime = isfield(data, 'time');
hasfreq = isfield(data, 'freq');

assert((hastime && hasfreq), 'please use ft_singleplotER for time-only or frequency-only data');

xparam = 'time';
yparam = 'freq';

% check whether rpt/subj is present and remove if necessary
dimord = getdimord(data, cfg.parameter);
dimtok = tokenize(dimord, '_');
hasrpt = any(ismember(dimtok, {'rpt' 'subj'}));

if ~hasrpt
  assert(isequal(cfg.trials, 'all') || isequal(cfg.trials, 1), 'incorrect specification of cfg.trials for data without repetitions');
else
  assert(~isempty(cfg.trials), 'empty specification of cfg.trials for data with repetitions');
end

% parse cfg.channel 
if isfield(cfg, 'channel') && isfield(data, 'label')
  cfg.channel = ft_channelselection(cfg.channel, data.label);
elseif isfield(cfg, 'channel') && isfield(data, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
end

% Apply baseline correction:
if ~strcmp(cfg.baseline, 'no')
  % keep mask-parameter if it is set
  if ~isempty(cfg.maskparameter)
    tempmask = data.(cfg.maskparameter);
  end
  data = ft_freqbaseline(cfg, data);
  % put mask-parameter back if it is set
  if ~isempty(cfg.maskparameter)
    data.(cfg.maskparameter) = tempmask;
  end
end

% channels should NOT be selected and averaged here, since a topoplot might follow in interactive mode
tmpcfg = keepfields(cfg, {'showcallinfo', 'trials'});
if hasrpt
  tmpcfg.avgoverrpt = 'yes';
else
  tmpcfg.avgoverrpt = 'no';
end
tmpvar = data;
[data] = ft_selectdata(tmpcfg, data);
% restore the provenance information and put back cfg.channel
tmpchannel  = cfg.channel;
[cfg, data] = rollback_provenance(cfg, data);
cfg.channel = tmpchannel;


if isfield(tmpvar, cfg.maskparameter) && ~isfield(data, cfg.maskparameter)
  % the mask parameter is not present after ft_selectdata, because it is
  % not included in all input arguments. Make the same selection and copy
  % it over
  tmpvar = ft_selectdata(tmpcfg, tmpvar);
  data.(cfg.maskparameter) = tmpvar.(cfg.maskparameter);
end

clear tmpvar tmpcfg dimord dimtok hastime hasfreq hasrpt

% ensure that the preproc specific options are located in the cfg.preproc
% substructure, but also ensure that the field 'refchannel' remains at the
% highest level in the structure. This is a little hack by JM because the field
% refchannel can relate to connectivity or to an EEg reference.

if isfield(cfg, 'refchannel'), refchannelincfg = cfg.refchannel; cfg = rmfield(cfg, 'refchannel'); end
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});
if exist('refchannelincfg', 'var'), cfg.refchannel  = refchannelincfg; end

if ~isempty(cfg.preproc)
  % preprocess the data, i.e. apply filtering, baselinecorrection, etc.
  fprintf('applying preprocessing options\n');
  if ~isfield(cfg.preproc, 'feedback')
    cfg.preproc.feedback = cfg.interactive;
  end
  data = ft_preprocessing(cfg.preproc, data);
end

% Handle the bivariate case
dimord = getdimord(data, cfg.parameter);
if startsWith(dimord, 'chan_chan_') || startsWith(dimord, 'chancmb_')
  % convert the bivariate data to univariate and call this plotting function again
  cfg.originalfunction = 'ft_singleplotTFR';
  cfg.trials = 'all'; % trial selection has been taken care off
  bivariate_common(cfg, data);
  return
end

% Apply channel-type specific scaling
tmpcfg = keepfields(cfg, {'parameter', 'chanscale', 'ecgscale', 'eegscale', 'emgscale', 'eogscale', 'gradscale', 'magscale', 'megscale', 'mychan', 'mychanscale'});
[data] = chanscale_common(tmpcfg, data);


%% Section 3: select the data to be plotted and determine min/max range

% Take the subselection of channels that is contained in the layout, this is the same in all datasets
[selchan] = match_str(data.label, cfg.channel);

% Get physical min/max range of x, i.e. time
if strcmp(cfg.xlim, 'maxmin')
  xmin = min(data.(xparam));
  xmax = max(data.(xparam));
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Get the index of the nearest bin
xminindx = nearest(data.(xparam), xmin);
xmaxindx = nearest(data.(xparam), xmax);
xmin = data.(xparam)(xminindx);
xmax = data.(xparam)(xmaxindx);
selx = xminindx:xmaxindx;
xval = data.(xparam)(selx);

% Get physical min/max range of y, i.e. frequency
if strcmp(cfg.ylim, 'maxmin')
  ymin = min(data.(yparam));
  ymax = max(data.(yparam));
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Get the index of the nearest bin
yminindx = nearest(data.(yparam), ymin);
ymaxindx = nearest(data.(yparam), ymax);
ymin = data.(yparam)(yminindx);
ymax = data.(yparam)(ymaxindx);
sely = yminindx:ymaxindx;
yval = data.(yparam)(sely);

% test if X and Y are linearly spaced (to within 10^-12): % FROM UIMAGE
dx = min(diff(xval));  % smallest interval for X
dy = min(diff(yval));  % smallest interval for Y
evenx = all(abs(diff(xval)/dx-1)<1e-12);     % true if X is linearly spaced
eveny = all(abs(diff(yval)/dy-1)<1e-12);     % true if Y is linearly spaced

if ~evenx || ~eveny
  ft_warning('(one of the) axis is/are not evenly spaced, but plots are made as if axis are linear')
end

% masking is only possible for evenly spaced axis
if strcmp(cfg.masknans, 'yes') && (~evenx || ~eveny)
  ft_warning('(one of the) axis are not evenly spaced -> nans cannot be masked out -> cfg.masknans is set to ''no'';')
  cfg.masknans = 'no';
end

% the usual data is chan_freq_time, but other dimords should also work
dimtok = tokenize(dimord, '_');
datamatrix = data.(cfg.parameter);
[c, ia, ib] = intersect(dimtok, {'chan', yparam, xparam});
datamatrix = permute(datamatrix, ia);
datamatrix = datamatrix(selchan, sely, selx);

if ~isempty(cfg.maskparameter)
  maskmatrix = data.(cfg.maskparameter)(selchan, sely, selx);
  if cfg.maskalpha ~= 1
    maskmatrix( maskmatrix) = 1;
    maskmatrix(~maskmatrix) = cfg.maskalpha;
  end
else
  % create an Nx0x0 matrix
  maskmatrix = zeros(length(selchan), 0, 0);
end

%% Section 4: do the actual plotting

cla
hold on

zval = mean(datamatrix, 1); % over channels
zval = reshape(zval, size(zval,2), size(zval,3));
mask = squeeze(mean(maskmatrix, 1)); % over channels

% Get physical z-axis range (color axis):
if strcmp(cfg.zlim, 'maxmin')
  zmin = nanmin(zval(:));
  zmax = nanmax(zval(:));
elseif strcmp(cfg.zlim, 'maxabs')
  zmin = -nanmax(abs(zval(:)));
  zmax =  nanmax(abs(zval(:)));
elseif strcmp(cfg.zlim, 'zeromax')
  zmin = 0;
  zmax = nanmax(zval(:));
elseif strcmp(cfg.zlim, 'minzero')
  zmin = nanmin(zval(:));
  zmax = 0;
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% Draw the data and mask NaN's if requested
if isequal(cfg.masknans, 'yes') && isempty(cfg.maskparameter)
  nans_mask = ~isnan(zval);
  mask = double(nans_mask);
  ft_plot_matrix(xval, yval, zval, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask)
elseif isequal(cfg.masknans, 'yes') && ~isempty(cfg.maskparameter)
  nans_mask = ~isnan(zval);
  mask = mask .* nans_mask;
  mask = double(mask);
  ft_plot_matrix(xval, yval, zval, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask)
elseif isequal(cfg.masknans, 'no') && ~isempty(cfg.maskparameter)
  mask = double(mask);
  ft_plot_matrix(xval, yval, zval, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask)
else
  ft_plot_matrix(xval, yval, zval, 'clim', [zmin zmax], 'tag', 'cip')
end

% set colormap
if isfield(cfg, 'colormap')
  if ~isnumeric(cfg.colormap)
    cfg.colormap = colormap(cfg.colormap);
  end
  if size(cfg.colormap,2)~=3
    ft_error('colormap must be a Nx3 matrix');
  else
    set(gcf, 'colormap', cfg.colormap);
  end
end

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

axis xy

if isequal(cfg.colorbar, 'yes')
  % tag the colorbar so we know which axes are colorbars
  colorbar('tag', 'ft-colorbar');
end

% Set callback to adjust color axis
if strcmp('yes', cfg.hotkeys)
  %  Attach data and cfg to figure and attach a key listener to the figure
  set(gcf, 'KeyPressFcn', {@key_sub, xmin, xmax, ymin, ymax, zmin, zmax})
end

% Create axis title containing channel name(s) and channel number(s):
if ~isempty(cfg.title)
  t = cfg.title;
else
  if length(cfg.channel) == 1
    t = [char(cfg.channel) ' / ' num2str(selchan) ];
  else
    t = sprintf('mean(%0s)', join_str(', ', cfg.channel));
  end
end
title(t, 'fontsize', cfg.fontsize);

% set the figure window title, add channel labels if number is small
if isempty(get(gcf, 'Name'))
  if length(selchan) < 5
    chans = join_str(', ', cfg.channel);
  else
    chans = '<multiple channels>';
  end
  if isempty(cfg.figurename)
    set(gcf, 'Name', sprintf('%d: %s: %s (%s)', double(gcf), mfilename, dataname, chans));
    set(gcf, 'NumberTitle', 'off');
  else
    set(gcf, 'name', cfg.figurename);
    set(gcf, 'NumberTitle', 'off');
  end
end

axis tight
hold off

% Make the figure interactive
if strcmp(cfg.interactive, 'yes')
  % add the cfg/data information to the figure under identifier linked to this axis
  ident             = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
  set(gca, 'tag',ident);
  info                  = guidata(gcf);
  info.(ident).dataname = dataname;
  info.(ident).cfg      = cfg;
  info.(ident).data     = data;
  guidata(gcf, info);
  set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR}, 'event', 'WindowButtonMotionFcn'});
end

% add a menu to the figure, but only if the current figure does not have subplots
% also, delete any possibly existing previous menu, this is safe because delete([]) does nothing
delete(findobj(gcf, 'type', 'uimenu', 'label', 'FieldTrip'));
if numel(findobj(gcf, 'type', 'axes', '-not', 'tag', 'ft-colorbar')) <= 1
  ftmenu = uimenu(gcf, 'Label', 'FieldTrip');
  uimenu(ftmenu, 'Label', 'Show pipeline',  'Callback', {@menu_pipeline, cfg});
  uimenu(ftmenu, 'Label', 'About',  'Callback', @menu_about);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance

if ~nargout
  % don't return anything
  clear cfg
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotTFR(range, varargin)
% fetch cfg/data based on axis indentifier given as tag
ident  = get(gca, 'tag');
info   = guidata(gcf);
cfg    = info.(ident).cfg;
data   = info.(ident).data;
if ~isempty(range)
  cfg = removefields(cfg, 'inputfile');   % the reading has already been done and varargin contains the data
  cfg = removefields(cfg, 'showlabels');  % this is not allowed in topoplotER
  cfg.trials = 'all';                     % trial selection has already been taken care of
  cfg.baseline = 'no';                    % make sure the next function does not apply a baseline correction again
  cfg.channel = 'all';                    % make sure the topo displays all channels, not just the ones in this singleplot
  cfg.comment = 'auto';
  cfg.dataname = info.(ident).dataname;   % put data name in here, this cannot be resolved by other means
  cfg.xlim = range(1:2);
  cfg.ylim = range(3:4);
  fprintf('selected cfg.xlim = [%f %f]\n', cfg.xlim(1), cfg.xlim(2));
  fprintf('selected cfg.ylim = [%f %f]\n', cfg.ylim(1), cfg.ylim(2));
  % ensure that the new figure appears at the same position
  f = figure('Position', get(gcf, 'Position'));
  ft_topoplotTFR(cfg, data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata, varargin)
xlimits = xlim;
ylimits = ylim;
climits = caxis;
incr_x = abs(xlimits(2) - xlimits(1)) /10;
incr_y = abs(ylimits(2) - ylimits(1)) /10;
incr_c = abs(climits(2) - climits(1)) /10;

if length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:}, 'control')
  % TRANSLATE by 10%
  switch eventdata.Key
    case 'pageup'
      caxis([min(caxis)+incr_c max(caxis)+incr_c]);
    case 'pagedown'
      caxis([min(caxis)-incr_c max(caxis)-incr_c]);
    case 'leftarrow'
      xlim([xlimits(1)+incr_x xlimits(2)+incr_x])
    case 'rightarrow'
      xlim([xlimits(1)-incr_x xlimits(2)-incr_x])
    case 'uparrow'
      ylim([ylimits(1)-incr_y ylimits(2)-incr_y])
    case 'downarrow'
      ylim([ylimits(1)+incr_y ylimits(2)+incr_y])
  end % switch
else
  % ZOOM by 10%
  switch eventdata.Key
    case 'pageup'
      caxis([min(caxis)-incr_c max(caxis)+incr_c]);
    case 'pagedown'
      caxis([min(caxis)+incr_c max(caxis)-incr_c]);
    case 'leftarrow'
      xlim([xlimits(1)-incr_x xlimits(2)+incr_x])
    case 'rightarrow'
      xlim([xlimits(1)+incr_x xlimits(2)-incr_x])
    case 'uparrow'
      ylim([ylimits(1)-incr_y ylimits(2)+incr_y])
    case 'downarrow'
      ylim([ylimits(1)+incr_y ylimits(2)-incr_y])
    case 'm'
      xlim([varargin{1} varargin{2}])
      ylim([varargin{3} varargin{4}])
      caxis([varargin{5} varargin{6}]);
  end % switch
end % if
