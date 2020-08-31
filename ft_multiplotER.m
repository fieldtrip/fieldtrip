function [cfg] = ft_multiplotER(cfg, varargin)

% FT_MULTIPLOTER plots the event-related potentials or event-related fields verus
% time, or the oscillatory activity (power or coherence) versus frequency. Multiple
% datasets can be overlayed. The plots are arranged according to their location
% specified in the layout.
%
% Use as
%   ft_multiplotER(cfg, data)
% or
%   ft_multiplotER(cfg, data, data2, ..., dataN)
%
% The data can be an event-related potential or field produced by
% FT_TIMELOCKANALYSIS, a power spectrum produced by FT_FREQANALYSIS or a coherence
% spectrum produced by FT_FREQDESCRIPTIVES.
%
% If you specify multiple datasets they should contain the same channels, etc.
%
% The configuration can have the following parameters:
%   cfg.parameter     = field to be plotted on y-axis, for example 'avg', 'powspctrm' or 'cohspctrm' (default is automatic)
%   cfg.maskparameter = field in the first dataset to be used for marking significant data
%   cfg.maskstyle     = style used for masking of data, 'box', 'thickness' or 'saturation' (default = 'box')
%   cfg.maskfacealpha = mask transparency value between 0 and 1
%   cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim          = 'maxmin', 'maxabs', 'zeromax', 'minzero', or [ymin ymax] (default = 'maxmin')
%   cfg.gradscale     = number, scaling to apply to the MEG gradiometer channels prior to display
%   cfg.magscale      = number, scaling to apply to the MEG magnetometer channels prior to display
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.refchannel    = name of reference channel for visualising connectivity, can be 'gui'
%   cfg.baseline      = 'yes', 'no' or [time1 time2] (default = 'no'), see FT_TIMELOCKBASELINE or FT_FREQBASELINE
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.axes          = string, 'yes' or 'no' whether to draw x- and y-axes for each graph (default = 'yes')
%   cfg.box           = string, 'yes' or 'no' whether to draw a box around each graph (default = 'no')
%   cfg.showlabels    = 'yes' or 'no' (default = 'no')
%   cfg.showoutline   = 'yes' or 'no' (default = 'no')
%   cfg.showscale     = 'yes' or 'no' (default = 'yes')
%   cfg.showcomment   = 'yes' or 'no' (default = 'yes')
%   cfg.comment       = string of text (default = date + limits)
%                       Add 'comment' to graph (according to COMNT in the layout)
%   cfg.limittext     = add user-defined text instead of cfg.comment, (default = cfg.comment)
%   cfg.fontsize      = font size of comment and labels (default = 8)
%   cfg.interactive   = 'yes' or 'no', make the plot interactive (default = 'yes')
%                       In a interactive plot you can select areas and produce a new
%                       interactive plot when a selected area is clicked. Multiple areas
%                       can be selected by holding down the SHIFT key.
%   cfg.renderer      = 'painters', 'zbuffer', ' opengl' or 'none' (default = [])
%   cfg.colorgroups   = 'sequential', 'allblack', 'labelcharN' (N = Nth character in label), 'chantype' or a vector
%                       with the length of the number of channels defining the groups (default = 'sequential')
%   cfg.linestyle     = linestyle/marker type, see options of the PLOT function (default = '-')
%                       can be a single style for all datasets, or a cell-array containing one style for each dataset
%   cfg.linewidth     = linewidth in points (default = 0.5)
%   cfg.linecolor     = color(s) used for plotting the dataset(s) (default = 'brgkywrgbkywrgbkywrgbkyw')
%                       alternatively, colors can be specified as Nx3 matrix of RGB values
%   cfg.directionality = '', 'inflow' or 'outflow' specifies for connectivity measures whether the
%                       inflow into a node, or the outflow from a node is plotted. The (default) behavior
%                       of this option depends on the dimord of the input data (see below).
%   cfg.layout        = specify the channel layout for plotting using one of the supported ways (see below).

% The following options for the scaling of the EEG, EOG, ECG, EMG, MEG and NIRS channels
% is optional and can be used to bring the absolute numbers of the different
% channel types in the same range (e.g. fT and uV). The channel types are determined
% from the input data using FT_CHANNELSELECTION.
%   cfg.eegscale      = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale      = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale      = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale      = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale      = number, scaling to apply to the MEG channels prior to display
%   cfg.gradscale     = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale      = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.nirsscale     = number, scaling to apply to the NIRS channels prior to display
%   cfg.mychanscale   = number, scaling to apply to the channels specified in cfg.mychan
%   cfg.mychan        = Nx1 cell-array with selection of channels
%   cfg.chanscale     = Nx1 vector with scaling factors, one per channel specified in cfg.channel
%
% For the plotting of directional connectivity data the cfg.directionality option
% determines what is plotted. The default value and the supported functionality
% depend on the dimord of the input data. If the input data is of dimord
% 'chan_chan_XXX', the value of directionality determines whether, given the
% reference channel(s), the columns (inflow), or rows (outflow) are selected for
% plotting. In this situation the default is 'inflow'. Note that for undirected
% measures, inflow and outflow should give the same output. If the input data is of
% dimord 'chancmb_XXX', the value of directionality determines whether the rows in
% data.labelcmb are selected. With 'inflow' the rows are selected if the
% refchannel(s) occur in the right column, with 'outflow' the rows are selected if
% the refchannel(s) occur in the left column of the labelcmb-field. Default in this
% case is '', which means that all rows are selected in which the refchannel(s)
% occur. This is to robustly support linearly indexed undirected connectivity
% metrics. In the situation where undirected connectivity measures are linearly
% indexed, specifying 'inflow' or 'outflow' can result in unexpected behavior.
%
% The layout defines how the channels are arranged and what the size of each
% subplot is. You can specify the layout in a variety of ways:
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
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure. For this particular function, the
% data should be provided as a cell-array.
%
% See also FT_MULTIPLOTTFR, FT_SINGLEPLOTER, FT_SINGLEPLOTTFR, FT_TOPOPLOTER,
% FT_TOPOPLOTTFR, FT_PREPARE_LAYOUT

% Undocumented local options:
% cfg.layoutname
% cfg.preproc
% cfg.orient = landscape/portrait

% Copyright (C) 2003-2006, Ole Jensen
% Copyright (C) 2007-2011, Roemer van der Meij & Jan-Mathijs Schoffelen
% Copyright (C) 2012-2017, F.C. Donders Centre
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

Ndata = length(varargin);
for i=1:Ndata
  % check if the input data is valid for this function
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'timelock', 'freq'});
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'cohrefchannel', 'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed', {'hlim', 'xlim'});
cfg = ft_checkconfig(cfg, 'renamed', {'matrixside',  'directionality'});
cfg = ft_checkconfig(cfg, 'renamed', {'vlim', 'ylim'});
cfg = ft_checkconfig(cfg, 'renamed', {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'renamed', {'graphcolor', 'linecolor'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedback', 'inflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim', 'absmax', 'maxabs'});
cfg = ft_checkconfig(cfg, 'unused',  {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamed', {'newfigure', 'figure'});
% cfg = ft_checkconfig(cfg, 'deprecated', {'xparam'});

% set the defaults
cfg.baseline       = ft_getopt(cfg, 'baseline', 'no');
cfg.trials         = ft_getopt(cfg, 'trials', 'all', 1);
cfg.xlim           = ft_getopt(cfg, 'xlim', 'maxmin');
cfg.ylim           = ft_getopt(cfg, 'ylim', 'maxmin');
cfg.comment        = ft_getopt(cfg, 'comment', date);
cfg.limittext      = ft_getopt(cfg, 'limittext', 'default');
cfg.axes           = ft_getopt(cfg, 'axes', 'yes');
cfg.showlabels     = ft_getopt(cfg, 'showlabels', 'no');
cfg.showoutline    = ft_getopt(cfg, 'showoutline', 'no');
cfg.showscale      = ft_getopt(cfg, 'showscale',   'yes');
cfg.showcomment    = ft_getopt(cfg, 'showcomment', 'yes');
cfg.box            = ft_getopt(cfg, 'box', 'no');
cfg.fontsize       = ft_getopt(cfg, 'fontsize', 8);
cfg.fontweight     = ft_getopt(cfg, 'fontweight');
cfg.interpreter    = ft_getopt(cfg, 'interpreter', 'none');  % none, tex or latex
cfg.interactive    = ft_getopt(cfg, 'interactive', 'yes');
cfg.orient         = ft_getopt(cfg, 'orient', 'landscape');
cfg.maskparameter  = ft_getopt(cfg, 'maskparameter');
cfg.colorgroups    = ft_getopt(cfg, 'colorgroups', 'condition');
cfg.linecolor      = ft_getopt(cfg, 'linecolor', 'brgkywrgbkywrgbkywrgbkyw');
cfg.linestyle      = ft_getopt(cfg, 'linestyle', '-');
cfg.linewidth      = ft_getopt(cfg, 'linewidth', 0.5);
cfg.maskstyle      = ft_getopt(cfg, 'maskstyle', 'box');
cfg.maskfacealpha  = ft_getopt(cfg, 'maskfacealpha', 1);
cfg.channel        = ft_getopt(cfg, 'channel', 'all');
cfg.directionality = ft_getopt(cfg, 'directionality', '');
cfg.figurename     = ft_getopt(cfg, 'figurename');
cfg.preproc        = ft_getopt(cfg, 'preproc');
cfg.frequency      = ft_getopt(cfg, 'frequency', 'all'); % needed for frequency selection with TFR data
cfg.latency        = ft_getopt(cfg, 'latency', 'all'); % needed for latency selection with TFR data, FIXME, probably not used
cfg.renderer       = ft_getopt(cfg, 'renderer'); % let MATLAB decide on the default

% check for linestyle being a cell-array
if ischar(cfg.linestyle)
  cfg.linestyle = repmat({cfg.linestyle}, 1, Ndata);
end
% check it's length, and lengthen it if does not have enough styles in it
if (length(cfg.linestyle) > 1) && (length(cfg.linestyle) < Ndata )
  ft_error('either specify cfg.linestyle as a cell-array with one cell for each dataset, or only specify one linestyle')
elseif (length(cfg.linestyle) == 1)
  cfg.linestyle = repmat(cfg.linestyle, 1, Ndata);
end

% this is needed for the figure title and correct labeling of linecolor later on
if isfield(cfg, 'dataname') && ~isempty(cfg.dataname)
  dataname = cfg.dataname;
elseif isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  dataname = cfg.inputfile;
elseif nargin>1
  dataname = arrayfun(@inputname, 2:nargin, 'UniformOutput', false);
else
  dataname = {};
end

%% Section 2: data handling, this also includes converting bivariate (chan_chan and chancmb) into univariate data

for i=1:Ndata
  dtype{i}   = ft_datatype(varargin{i});
  hastime(i) = isfield(varargin{i}, 'time');
  hasfreq(i) = isfield(varargin{i}, 'freq');
end

% check if the input has consistent datatypes
if ~all(strcmp(dtype, dtype{1})) || ~all(hastime==hastime(1)) || ~all(hasfreq==hasfreq(1))
  ft_error('different datatypes are not allowed as input');
end
dtype   = dtype{1};
hastime = hastime(1);
hasfreq = hasfreq(1);

% Set x/y/parameter defaults according to datatype and dimord
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
      xparam = 'time';  % further down the code will compute the average over selected frequencies
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
  tmpcfg = keepfields(cfg, {'baseline', 'baselinetype', 'baselinewindow', 'demean', 'parameter', 'channel'});
  for i=1:Ndata
    % keep mask-parameter if it is set
    if ~isempty(cfg.maskparameter)
      tempmask = varargin{i}.(cfg.maskparameter);
    end
    if strcmp(dtype, 'timelock') && strcmp(xparam, 'time')
      varargin{i} = ft_timelockbaseline(tmpcfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(xparam, 'time')
      varargin{i} = ft_freqbaseline(tmpcfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(xparam, 'freq')
      ft_error('Baseline correction is not supported for spectra without a time dimension');
    else
      ft_warning('Baseline correction not applied, please set xparam');
    end
    % put mask-parameter back if it is set
    if ~isempty(cfg.maskparameter)
      varargin{i}.(cfg.maskparameter) = tempmask;
    end
  end
end

% channels SHOULD be selected here, as no interactive action produces a new multiplot
tmpcfg = keepfields(cfg, {'channel', 'showcallinfo', 'trials'});
if hasrpt
  tmpcfg.avgoverrpt = 'yes';
else
  tmpcfg.avgoverrpt = 'no';
end
if hastime && hasfreq
  tmpcfg.avgoverfreq = 'yes';       % average over the selected frequencies
  tmpcfg.frequency = cfg.frequency; % not to be confused with cfg.xlim or cfg.ylim
  tmpcfg.keepfreqdim = 'no';
else
  tmpcfg.avgoverfreq = 'no';
end
tmpvar = varargin{1};
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

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

% handle the bivariate case
dimord = getdimord(varargin{1}, cfg.parameter);
if startsWith(dimord, 'chan_chan_') || startsWith(dimord, 'chancmb_')
  % convert the bivariate data to univariate and call this plotting function with univariate input
  cfg.originalfunction = 'ft_multiplotER';
  cfg.trials = 'all'; % trial selection has been taken care off
  bivariate_common(cfg, varargin{:});
  return
end

% Apply channel-type specific scaling
fn = fieldnames(cfg);
fn = setdiff(fn, {'skipscale', 'showscale', 'gridscale'}); % these are for the layout and plotting, not for CHANSCALE_COMMON
fn = fn(endsWith(fn, 'scale') | startsWith(fn, 'mychan') | strcmp(fn, 'channel') | strcmp(fn, 'parameter'));
tmpcfg = keepfields(cfg, fn);
if ~isempty(tmpcfg)
  for i=1:Ndata
    varargin{i} = chanscale_common(tmpcfg, varargin{i});
  end
  % remove the scaling fields from the configuration, to prevent them from being called again in interactive mode
  % but keep the parameter and channel field
  cfg = removefields(cfg, setdiff(fn, {'parameter', 'channel'}));
else
  % do nothing
end

%% Section 3: select the data to be plotted and determine min/max range

% Read or create the layout that will be used for plotting
tmpcfg = keepfields(cfg, {'layout', 'rows', 'columns', 'commentpos', 'skipcomnt', 'scalepos', 'skipscale', 'projection', 'viewpoint', 'rotate', 'width', 'height', 'elec', 'grad', 'opto', 'showcallinfo'});
cfg.layout = ft_prepare_layout(tmpcfg, varargin{1});

% Take the subselection of channels that is contained in the layout, this is the same in all datasets
[selchan, sellay] = match_str(varargin{1}.label, cfg.layout.label);

% Get physical min/max range of x, i.e. time or frequency
if strcmp(cfg.xlim, 'maxmin')
  % Find maxmin throughout all varargins:
  xmin = [];
  xmax = [];
  for i=1:Ndata
    xmin = nanmin([xmin varargin{i}.(xparam)]);
    xmax = nanmax([xmax varargin{i}.(xparam)]);
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

% Get physical y-axis range, i.e. of the parameter to be plotted
if ~isnumeric(cfg.ylim)
  % Find maxmin throughout all varargins
  ymin = +inf;
  ymax = -inf;
  for i=1:Ndata
    % Select the channels in the data that match with the layout and that are selected for plotting
    dat = varargin{i}.(cfg.parameter)(selchan,selx);
    ymin = min(ymin, min(dat(:)));
    ymax = max(ymax, max(dat(:)));
  end
  switch cfg.ylim
    case 'maxmin'
      % keep them as they are
    case 'maxabs'
      ymax = max(abs(ymax), abs(ymin));
      ymin = -ymax;
    case 'zeromax'
      ymin = 0;
    case 'minzero'
      ymax = 0;
    otherwise
      ft_error('invalid specification of cfg.ylim');
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
  % one value for each channel-time point
  maskmatrix = varargin{1}.(cfg.maskparameter)(selchan, selx);
else
  % create an Nx0 matrix
  maskmatrix = zeros(length(selchan), 0);
end

chanX      = cfg.layout.pos(sellay, 1);
chanY      = cfg.layout.pos(sellay, 2);
chanWidth  = cfg.layout.width(sellay);
chanHeight = cfg.layout.height(sellay);
chanLabel  = cfg.layout.label(sellay);

%% Section 4: do the actual plotting

% determine the coloring of channels/conditions
linecolor = linecolor_common(cfg, varargin{:});

% open a new figure, or add it to the existing one
open_figure(keepfields(cfg, {'figure', 'clearfigure', 'position', 'visible', 'renderer', 'figurename', 'title'}));

if ischar(linecolor)
  set(gca, 'ColorOrder', char2rgb(linecolor))
elseif isnumeric(linecolor)
  set(gca, 'ColorOrder', linecolor)
end

% Plot the data
for m=1:length(selchan)
  mask = maskmatrix(m, :);
  if strcmp(cfg.maskstyle, 'difference')
    % combine the conditions in a single plot, highlight the difference
    yval = squeeze(datamatrix(:,m,:));
    % Clip out of bounds y values:
    yval(yval > ymax) = ymax;
    yval(yval < ymin) = ymin;
    ft_plot_vector(xval, yval, 'width', chanWidth(m), 'height', chanHeight(m), 'hpos', chanX(m), 'vpos', chanY(m), 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'color', linecolor, 'style', cfg.linestyle{1}, 'linewidth', cfg.linewidth, 'axis', cfg.axes, 'highlight', mask, 'highlightstyle', cfg.maskstyle, 'facealpha', cfg.maskfacealpha);
  else
    % loop over the conditions, plot them on top of each other
    for i=1:Ndata
      yval = squeeze(datamatrix(i,m,:));
      % clip out of bounds y values:
      yval(yval > ymax) = ymax;
      yval(yval < ymin) = ymin;
      % select the color for the channel/condition
      if strcmp(cfg.colorgroups, 'condition')
        color = linecolor(i,:);
      else
        color = linecolor(m,:);
      end
      ft_plot_vector(xval, yval, 'width', chanWidth(m), 'height', chanHeight(m), 'hpos', chanX(m), 'vpos', chanY(m), 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'color', color, 'style', cfg.linestyle{i}, 'linewidth', cfg.linewidth, 'axis', cfg.axes, 'highlight', mask, 'highlightstyle', cfg.maskstyle, 'facealpha', cfg.maskfacealpha);
    end
  end
end % for number of channels

% plot the layout, labels and outline
ft_plot_layout(cfg.layout, 'box', cfg.box, 'label', cfg.showlabels, 'outline', cfg.showoutline, 'point', 'no', 'mask', 'no', 'fontsize', cfg.fontsize, 'labelyoffset', 1.4*median(cfg.layout.height/2), 'labelalignh', 'center', 'chanindx', find(~ismember(cfg.layout.label, {'COMNT', 'SCALE'})), 'interpreter', cfg.interpreter);

% write comment
if istrue(cfg.showcomment)
  % Add the colors of the different conditions/datasets to the comment
  colorLabels = [];
  if Ndata > 1 && strcmp(cfg.colorgroups, 'condition')
    for i=1:Ndata
      if ischar(linecolor)
        colorLabels = [colorLabels '\n' dataname{i} '='         linecolor(i)     ];
      elseif isnumeric(linecolor)
        colorLabels = [colorLabels '\n' dataname{i} '=' num2str(linecolor(i, :)) ];
      end
    end
  end
  cfg.comment = [cfg.comment colorLabels];
  
  k = find(strcmp('COMNT', cfg.layout.label));
  if ~isempty(k)
    limittext = cfg.limittext;
    if ~strcmp(limittext, 'default')
      comment = limittext;
    else
      comment = cfg.comment;
      comment = sprintf('%0s\nxlim=[%.3g %.3g]', comment, xmin, xmax);
      comment = sprintf('%0s\nylim=[%.3g %.3g]', comment, ymin, ymax);
    end
    ft_plot_text(cfg.layout.pos(k, 1), cfg.layout.pos(k, 2), sprintf(comment), 'FontSize', cfg.fontsize, 'FontWeight', cfg.fontweight);
    % plot an invisible box, the text itself is not sufficient to get the automatic scaling of the figures axes to include COMNT
    xy(1) = cfg.layout.pos(k, 1) - cfg.layout.width(k, 1)/2;
    xy(2) = cfg.layout.pos(k, 1) + cfg.layout.width(k, 1)/2;
    xy(3) = cfg.layout.pos(k, 2) - cfg.layout.height(k, 1)/2;
    xy(4) = cfg.layout.pos(k, 2) + cfg.layout.height(k, 1)/2;
    ft_plot_box(xy, 'edgecolor', 'none');
  end
end

% Plot scales
if istrue(cfg.showscale)
  l = find(strcmp(cfg.layout.label, 'SCALE'));
  if ~isempty(l)
    x = cfg.layout.pos(l,1);
    y = cfg.layout.pos(l,2);
    plotScales([xmin xmax], [ymin ymax], x, y, chanWidth(1), chanHeight(1), cfg)
  end
end

axis tight
axis off

% Make the axis a little wider when boxes are shown
if strcmp(cfg.box, 'yes')
  abc = axis;
  axis(abc + [-1 +1 -1 +1]*mean(abs(abc))/10)
end

% Set orientation for printing if specified
if ~isempty(cfg.orient)
  orient(gcf, cfg.orient);
end

% set the figure window title
if ~isempty(dataname)
  set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), mfilename, join_str(', ', dataname)));
else
  set(gcf, 'Name', sprintf('%d: %s', double(gcf), mfilename));
end
set(gcf, 'NumberTitle', 'off');

% Make the figure interactive
if strcmp(cfg.interactive, 'yes')
  % add the cfg/data/channel information to the figure under identifier linked to this axis
  ident                 = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
  set(gca, 'tag', ident);
  info                  = guidata(gcf);
  info.(ident).x        = cfg.layout.pos(:, 1);
  info.(ident).y        = cfg.layout.pos(:, 2);
  info.(ident).label    = cfg.layout.label;
  info.(ident).dataname = dataname;
  info.(ident).cfg      = cfg;
  info.(ident).varargin = varargin;
  guidata(gcf, info);
  
  set(gcf, 'WindowButtonUpFcn',  {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER}, 'event', 'WindowButtonMotionFcn'});
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance
ft_postamble savefig

% add a menu to the figure, but only if the current figure does not have subplots
menu_fieldtrip(gcf, cfg, false);

if ~ft_nargout
  % don't return anything
  clear cfg
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotScales(hlim, vlim, hpos, vpos, width, height, cfg)

% the placement of all elements is identical
placement = {'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim};

ft_plot_box([hlim vlim], placement{:}, 'edgecolor', 'k');

if hlim(1)<=0 && hlim(2)>=0
  ft_plot_line([0 0], vlim, placement{:}, 'color', 'k');
  ft_plot_text(0, vlim(1), '0  ', placement{:}, 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', cfg.fontsize);
end

if vlim(1)<=0 && vlim(2)>=0
  ft_plot_line(hlim, [0 0], placement{:}, 'color', 'k');
  ft_plot_text(hlim(1), 0, '0  ', placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'middle', 'FontSize', cfg.fontsize);
end

ft_plot_text(hlim(1), vlim(1), [num2str(hlim(1), 3) ' '], placement{:}, 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',    'FontSize', cfg.fontsize);
ft_plot_text(hlim(2), vlim(1), [num2str(hlim(2), 3) ' '], placement{:}, 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', cfg.fontsize);
ft_plot_text(hlim(1), vlim(1), [num2str(vlim(1), 3) ' '], placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom', 'FontSize', cfg.fontsize);
ft_plot_text(hlim(1), vlim(2), [num2str(vlim(2), 3) ' '], placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'top',    'FontSize', cfg.fontsize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label, varargin)
% fetch cfg/data based on axis indentifier given as tag
ident       = get(gca,'tag');
info        = guidata(gcf);
cfg         = info.(ident).cfg;
datvarargin = info.(ident).varargin;
if ~isempty(label)
  cfg = removefields(cfg, 'inputfile');   % the reading has already been done and varargin contains the data
  cfg.baseline = 'no';                    % make sure the next function does not apply a baseline correction again
  cfg.channel = label;
  cfg.dataname = info.(ident).dataname;   % put data name in here, this cannot be resolved by other means
  cfg.trials = 'all';                     % trial selection has already been taken care of
  fprintf('selected cfg.channel = {%s}\n', join_str(', ', cfg.channel));
  % ensure that the new figure appears at the same position
  f = figure('position', get(gcf, 'Position'));
  ft_singleplotER(cfg, datvarargin{:});
end
