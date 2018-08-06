function [cfg] = ft_singleplotER(cfg, varargin)

% FT_SINGLEPLOTER plots the event-related fields or potentials of a single
% channel or the average over multiple channels. Multiple datasets can be
% overlayed.
%
% Use as
%   ft_singleplotER(cfg, data)
% or
%   ft_singleplotER(cfg, data1, data2, ..., datan)
%
% The data can be an erp/erf produced by FT_TIMELOCKANALYSIS, a power
% spectrum produced by FT_FREQANALYSIS or connectivity spectrum produced by
% FT_CONNECTIVITYANALYSIS.
%
% The configuration can have the following parameters:
%   cfg.parameter     = field to be plotted on y-axis (default depends on data.dimord)
%                       'avg', 'powspctrm' or 'cohspctrm'
%   cfg.maskparameter = field in the first dataset to be used for masking of data
%                       (not possible for mean over multiple channels, or when input contains multiple subjects
%                       or trials)
%   cfg.maskstyle     = style used for masking of data, 'box', 'thickness' or 'saturation' (default = 'box')
%   cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim          = 'maxmin', 'maxabs', 'zeromax', 'minzero', or [ymin ymax] (default = 'maxmin')
%   cfg.channel       = nx1 cell-array with selection of channels (default = 'all')
%                       see ft_channelselection for details
%   cfg.title          = string, title of plot
%   cfg.refchannel    = name of reference channel for visualising connectivity, can be 'gui'
%   cfg.baseline      = 'yes', 'no' or [time1 time2] (default = 'no'), see ft_timelockbaseline
%   cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
%   cfg.trials        = 'all' or a selection given as a 1xn vector (default = 'all')
%   cfg.fontsize      = font size of title (default = 8)
%   cfg.hotkeys       = enables hotkeys (leftarrow/rightarrow/uparrow/downarrow/m) for dynamic zoom and translation (ctrl+) of the axes
%   cfg.interactive   = interactive plot 'yes' or 'no' (default = 'yes')
%                       in a interactive plot you can select areas and produce a new
%                       interactive plot when a selected area is clicked. multiple areas
%                       can be selected by holding down the shift key.
%   cfg.renderer      = 'painters', 'zbuffer', ' opengl' or 'none' (default = [])
%   cfg.linestyle     = linestyle/marker type, see options of the PLOT function (default = '-')
%                       can be a single style for all datasets, or a cell-array containing one style for each dataset
%   cfg.linewidth     = linewidth in points (default = 0.5)
%   cfg.graphcolor    = color(s) used for plotting the dataset(s) (default = 'brgkywrgbkywrgbkywrgbkyw')
%                       alternatively, colors can be specified as nx3 matrix of rgb values
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
% to facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% if you specify this option the input data will be read from a *.mat
% file on disk. this mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SINGLEPLOTTFR, FT_MULTIPLOTER, FT_MULTIPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR

% Undocumented local options:
% cfg.zlim/xparam (set to a specific frequency range or time range [zmax zmin] for an average
% over the frequency/time bins for TFR data.  Use in conjunction with e.g. xparam = 'time', and cfg.parameter = 'powspctrm').

% Copyright (C) 2003-2006, Ole Jensen
% Copyright (C) 2006-2017, F.C. Donders Centre
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
ft_preamble loadvar varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

Ndata = numel(varargin);
for i=1:Ndata
  % check if the input data is valid for this function
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'timelock', 'freq'});
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'unused',     {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim', 'absmax', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedback',    'inflow'});
cfg = ft_checkconfig(cfg, 'renamed',    {'matrixside',     'directionality'});
cfg = ft_checkconfig(cfg, 'renamed',    {'channelindex',   'channel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'channelname',    'channel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'cohrefchannel',  'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'zparam',         'parameter'});
cfg = ft_checkconfig(cfg, 'deprecated', {'xparam'});

% set the defaults
cfg.baseline        = ft_getopt(cfg, 'baseline',      'no');
cfg.trials          = ft_getopt(cfg, 'trials',        'all', 1);
cfg.xlim            = ft_getopt(cfg, 'xlim',          'maxmin');
cfg.ylim            = ft_getopt(cfg, 'ylim',          'maxmin');
cfg.zlim            = ft_getopt(cfg, 'zlim',          'maxmin');
cfg.comment         = ft_getopt(cfg, 'comment',        strcat([date '\n']));
cfg.axes            = ft_getopt(cfg, ' axes',         'yes');
cfg.fontsize        = ft_getopt(cfg, 'fontsize',       8);
cfg.graphcolor      = ft_getopt(cfg, 'graphcolor',    'brgkywrgbkywrgbkywrgbkyw');
cfg.hotkeys         = ft_getopt(cfg, 'hotkeys',       'yes');
cfg.interactive     = ft_getopt(cfg, 'interactive',   'yes');
cfg.renderer        = ft_getopt(cfg, 'renderer',       []);
cfg.maskparameter   = ft_getopt(cfg, 'maskparameter',  []);
cfg.linestyle       = ft_getopt(cfg, 'linestyle',     '-');
cfg.linewidth       = ft_getopt(cfg, 'linewidth',      0.5);
cfg.maskstyle       = ft_getopt(cfg, 'maskstyle',     'box');
cfg.channel         = ft_getopt(cfg, 'channel',       'all');
cfg.title           = ft_getopt(cfg, 'title',          []);
cfg.directionality  = ft_getopt(cfg, 'directionality', []);
cfg.figurename      = ft_getopt(cfg, 'figurename',     []);
cfg.preproc         = ft_getopt(cfg, 'preproc',        []);
cfg.frequency       = ft_getopt(cfg, 'frequency',     'all'); % needed for frequency selection with TFR data
cfg.latency         = ft_getopt(cfg, 'latency',       'all'); % needed for latency selection with TFR data, FIXME, probably not used

if ischar(cfg.graphcolor)
  graphcolor = cfg.graphcolor(1:Ndata);
elseif isnumeric(cfg.graphcolor)
  graphcolor = cfg.graphcolor(1:Ndata,:);
end

% check for linestyle being a cell-array
if ischar(cfg.linestyle)
  cfg.linestyle = {cfg.linestyle};
end

% check linestyle length, and lengthen it if does not have enough styles in it
if (length(cfg.linestyle) < Ndata ) && (length(cfg.linestyle) > 1)
  ft_error('either specify cfg.linestyle as a cell-array with one cell for each dataset, or only specify one linestyle')
elseif (length(cfg.linestyle) < Ndata ) && (length(cfg.linestyle) == 1)
  tmpstyle = cfg.linestyle{1};
  cfg.linestyle = cell(Ndata , 1);
  for idataset = 1:Ndata
    cfg.linestyle{idataset} = tmpstyle;
  end
end

% this is needed for the figure title and correct labeling of graphcolor later on
if nargin>1
  if isfield(cfg, 'dataname')
    dataname = cfg.dataname;
  else
    dataname = cell(1,Ndata);
    for i=1:Ndata
      if ~isempty(inputname(i+1))
        dataname{i} = inputname(i+1);
      else
        dataname{i} = ['data' num2str(i,'%02d')];
      end
    end
  end
else  % data provided through cfg.inputfile
  dataname = cfg.inputfile;
end

%% Section 2: data handling, this also includes converting bivariate (chan_chan and chancmb) into univariate data

for i=1:Ndata
  dtype{i}    = ft_datatype(varargin{i});
  hastime(i)  = isfield(varargin{i}, 'time');
  hasfreq(i)  = isfield(varargin{i}, 'freq');
end

% check if the input has consistent datatypes
if ~all(strcmp(dtype, dtype{1})) || ~all(hastime==hastime(1)) || ~all(hasfreq==hasfreq(1))
  ft_error('different datatypes are not allowed as input');
else
  dtype   = dtype{1};
  hastime = hastime(1);
  hasfreq = hasfreq(1);
end

% Set x/y/parameter according to datatype and dimord
switch dtype
  case 'timelock'
    xparam = 'time';
    if isfield(varargin{1}, 'trial')
      cfg.parameter = ft_getopt(cfg, 'parameter', 'trial');
    elseif isfield(varargin{1}, 'individual')
      cfg.parameter = ft_getopt(cfg, 'parameter', 'individual');
    elseif isfield(varargin{1}, 'avg')
      cfg.parameter = ft_getopt(cfg, 'parameter', 'avg');
    end
  case 'freq'
    if hastime && hasfreq
      xparam = 'time';   % average over selected frequencies
    else
      xparam = 'freq';
    end
    cfg.parameter = ft_getopt(cfg, 'parameter', 'powspctrm');
  case 'comp'
    % not supported
  otherwise
    % not supported
end

% check whether rpt/subj is present and remove if necessary
dimord = getdimord(varargin{1}, cfg.parameter);
dimtok = tokenize(dimord, '_');
hasrpt = any(ismember(dimtok, {'rpt' 'subj'}));

if ~hasrpt
  assert(isequal(cfg.trials, 'all') || isequal(cfg.trials, 1), 'incorrect specification of cfg.trials for data without repetitions');
else
  assert(~isempty(cfg.trials), 'empty specification of cfg.trials for data with repetitions');
end

% parse cfg.channel
if isfield(cfg, 'channel') && isfield(varargin{1}, 'label')
  cfg.channel = ft_channelselection(cfg.channel, varargin{1}.label);
elseif isfield(cfg, 'channel') && isfield(varargin{1}, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(varargin{1}.labelcmb(:)));
end

% apply baseline correction
if ~strcmp(cfg.baseline, 'no')
  for i=1:Ndata
    % keep mask-parameter if it is set
    if ~isempty(cfg.maskparameter)
      tempmask = varargin{i}.(cfg.maskparameter);
    end
    if strcmp(dtype, 'timelock') && strcmp(xparam, 'time')
      varargin{i} = ft_timelockbaseline(cfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(xparam, 'time')
      varargin{i} = ft_freqbaseline(cfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(xparam, 'freq')
      ft_error('baseline correction is not supported for spectra without a time dimension');
    else
      ft_warning('baseline correction not applied, please set xparam');
    end
    % put mask-parameter back if it is set
    if ~isempty(cfg.maskparameter)
      varargin{i}.(cfg.maskparameter) = tempmask;
    end
  end
end


% channels should NOT be selected and averaged here, since a topoplot might follow in interactive mode
tmpcfg = keepfields(cfg, {'showcallinfo', 'trials'});
if hasrpt
  tmpcfg.avgoverrpt = 'yes';
else
  tmpcfg.avgoverrpt = 'no';
end
if hastime && hasfreq
  tmpcfg.avgoverfreq = 'yes';       % average over selected frequencies
  tmpcfg.frequency = cfg.frequency; % not to be confused with cfg.xlim or cfg.ylim
  tmpcfg.keepfreqdim = 'no';
else
  tmpcfg.avgoverfreq = 'no';
end
tmpvar = varargin{1};
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information and put back cfg.channel
tmpchannel  = cfg.channel;
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
cfg.channel = tmpchannel;

if isfield(tmpvar, cfg.maskparameter) && ~isfield(varargin{1}, cfg.maskparameter)
  % the mask parameter is not present after ft_selectdata, because it is
  % not included in all input arguments. Make the same selection and copy
  % it over
  tmpvar = ft_selectdata(tmpcfg, tmpvar);
  varargin{1}.(cfg.maskparameter) = tmpvar.(cfg.maskparameter);
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
  for i=1:Ndata
    varargin{i} = ft_preprocessing(cfg.preproc, varargin{i});
  end
end

% Handle the bivariate case
dimord = getdimord(varargin{1}, cfg.parameter);
if startsWith(dimord, 'chan_chan_') || startsWith(dimord, 'chancmb_')
  % convert the bivariate data to univariate and call this plotting function again
  cfg.originalfunction = 'ft_singleplotER';
  cfg.trials = 'all'; % trial selection has been taken care off
  bivariate_common(cfg, varargin{:});
  return
end

% Apply channel-type specific scaling
tmpcfg = keepfields(cfg, {'parameter', 'chanscale', 'ecgscale', 'eegscale', 'emgscale', 'eogscale', 'gradscale', 'magscale', 'megscale', 'mychan', 'mychanscale'});
for i=1:Ndata
  varargin{i}= chanscale_common(tmpcfg, varargin{i});
end


%% Section 3: select the data to be plotted and determine min/max range

% Take the desided subselection of channels, this is the same in all datasets
[selchan] = match_str(varargin{1}.label, cfg.channel);

% Get physical min/max range of x, i.e. time or frequency
if strcmp(cfg.xlim, 'maxmin')
  % Find maxmin throughout all varargins:
  xmin = [];
  xmax = [];
  for i=1:Ndata
    xmin = min([xmin varargin{i}.(xparam)]);
    xmax = max([xmax varargin{i}.(xparam)]);
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Get the index of the nearest bin, this is the same in all datasets
xminindx = nearest(varargin{1}.(xparam), xmin);
xmaxindx = nearest(varargin{1}.(xparam), xmax);
xmin = varargin{1}.(xparam)(xminindx);
xmax = varargin{1}.(xparam)(xmaxindx);
selx = xminindx:xmaxindx;
xval = varargin{1}.(xparam)(selx);

% Get physical y-axis range, i.e. parameter to be plotted
if ~isnumeric(cfg.ylim)
  % Find maxmin throughout all varargins
  ymin = [];
  ymax = [];
  for i=1:Ndata
    % Select the channels in the data that match with the layout and that are selected for plotting
    dat = nanmean(varargin{i}.(cfg.parameter)(selchan,selx),1); % mean over channels, as that is what will be plotted
    ymin = nanmin([ymin nanmin(nanmin(nanmin(dat)))]);
    ymax = nanmax([ymax nanmax(nanmax(nanmax(dat)))]);
  end
  if strcmp(cfg.ylim, 'maxabs') % handle maxabs, make y-axis center on 0
    ymax = max([abs(ymax) abs(ymin)]);
    ymin = -ymax;
  elseif strcmp(cfg.ylim, 'zeromax')
    ymin = 0;
  elseif strcmp(cfg.ylim, 'minzero')
    ymax = 0;
  end
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Gather the data from all input data structures
datamatrix = zeros(Ndata, length(selchan), length(selx));
for i=1:Ndata
  datamatrix(i,:,:) = varargin{i}.(cfg.parameter)(selchan, selx);
end

if ~isempty(cfg.maskparameter)
  % one value for each channel, or one value for each channel-time point
  maskmatrix = varargin{1}.(cfg.maskparameter)(selchan, selx);
else
  % create an Nx0 matrix
  maskmatrix = zeros(length(selchan), 0);
end

%% Section 4: do the actual plotting

cla
hold on

yval = mean(datamatrix, 2); % over channels
yval = reshape(yval, size(yval,1), size(yval,3));
mask = squeeze(mean(maskmatrix, 1)); % over channels

ft_plot_vector(xval, yval, 'style', cfg.linestyle{i}, 'color', graphcolor, ...
  'highlight', mask, 'highlightstyle', cfg.maskstyle, 'linewidth', cfg.linewidth, ...
  'hlim', [xmin xmax], 'vlim', [ymin ymax]);

colorLabels = [];
if Ndata > 1
  for i=1:Ndata
    if ischar(graphcolor)
      colorLabels = [colorLabels '\n' dataname{i} '='         graphcolor(i)     ];
    elseif isnumeric(graphcolor)
      colorLabels = [colorLabels '\n' dataname{i} '=' num2str(graphcolor(i, :)) ];
    end
  end
end

% set xlim and ylim
if xmin~=xmax
  xlim([xmin xmax]);
end
if ymin~=ymax
  ylim([ymin ymax]);
end

% adjust mask box extents to ymin/ymax
if ~isempty(cfg.maskparameter)
  ptchs = findobj(gcf, 'type', 'patch');
  for i = 1:length(ptchs)
    YData = get(ptchs(i), 'YData');
    YData(YData == min(YData)) = ymin;
    YData(YData == max(YData)) = ymax;
    set(ptchs(i), 'YData',YData);
  end
end

% Set callback to adjust axes
if strcmp('yes', cfg.hotkeys)
  %  attach data and cfg to figure and attach a key listener to the figure
  set(gcf, 'KeyPressFcn', {@key_sub, xmin, xmax, ymin, ymax})
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
    if iscell(dataname)
      set(gcf, 'Name', sprintf('%d: %s: %s (%s)', double(gcf), mfilename, join_str(', ', dataname), chans));
    else
      set(gcf, 'Name', sprintf('%d: %s: %s (%s)', double(gcf), mfilename, dataname, chans));
    end
    set(gcf, 'NumberTitle', 'off');
  else
    set(gcf, 'name', cfg.figurename);
    set(gcf, 'NumberTitle', 'off');
  end
end

% set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

hold off

% Make the figure interactive
if strcmp(cfg.interactive, 'yes')
  % add the cfg/data/channel information to the figure under identifier linked to this axis
  ident                  = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
  set(gca, 'tag', ident);
  info                   = guidata(gcf);
  info.(ident).cfg       = cfg;
  info.(ident).varargin  = varargin;
  info.(ident).dataname  = dataname;
  guidata(gcf, info);
  set(gcf, 'windowbuttonupfcn',     {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER}, 'event', 'windowbuttonupfcn'});
  set(gcf, 'windowbuttondownfcn',   {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER}, 'event', 'windowbuttondownfcn'});
  set(gcf, 'windowbuttonmotionfcn', {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER}, 'event', 'windowbuttonmotionfcn'});
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
ft_postamble previous varargin
ft_postamble provenance

if ~nargout
  % don't return anything
  clear cfg
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotER(range, varargin)

% fetch cfg/data based on axis indentifier given as tag
ident    = get(gca, 'tag');
info     = guidata(gcf);
cfg      = info.(ident).cfg;
varargin = info.(ident).varargin;
if ~isempty(range)
  cfg = removefields(cfg, 'inputfile');   % the reading has already been done and varargin contains the data
  cfg = removefields(cfg, 'showlabels');  % this is not allowed in topoplotER
  cfg.baseline = 'no';                    % make sure the next function does not apply a baseline correction again
  cfg.dataname = info.(ident).dataname;   % put data name in here, this cannot be resolved by other means
  cfg.channel = 'all';                    % make sure the topo displays all channels, not just the ones in this singleplot
  cfg.comment = 'auto';
  cfg.trials = 'all';                     % trial selection has already been taken care of
  cfg.xlim = range(1:2);
  % if user specified a ylim, copy it over to the zlim of topoplot
  if isfield(cfg, 'ylim')
    cfg.zlim = cfg.ylim;
    cfg = rmfield(cfg, 'ylim');
  end
  fprintf('selected cfg.xlim = [%f %f]\n', cfg.xlim(1), cfg.xlim(2));
  % ensure that the new figure appears at the same position
  f = figure('position', get(gcf, 'Position'));
  ft_topoplotER(cfg, varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata, varargin)
xlimits = xlim;
ylimits = ylim;
incr_x = abs(xlimits(2) - xlimits(1)) /10;
incr_y = abs(ylimits(2) - ylimits(1)) /10;

if length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:}, 'control')
  % TRANSLATE by 10%
  switch eventdata.Key
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
  end % switch
end % if
