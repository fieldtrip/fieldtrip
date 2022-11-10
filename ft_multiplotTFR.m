function [cfg] = ft_multiplotTFR(cfg, data)

% FT_MULTIPLOTTFR plots the time-frequency representations of power or coherence
% in a topographical layout. The plots of the indivual sensors are arranged
% according to their location specified in the layout.
%
% Use as
%   ft_multiplotTFR(cfg, data)
%
% The data can be a time-frequency representation of power or coherence
% that was computed using the FT_FREQANALYSIS or FT_FREQDESCRIPTIVES
% functions.
%
% The configuration can have the following parameters:
%   cfg.parameter        = field to be represented as color (default depends on data.dimord)
%                          'powspctrm' or 'cohspctrm'
%   cfg.maskparameter    = field in the data to be used for masking of data, can be logical (e.g. significant data points) or numerical (e.g. t-values).
%                        (not possible for mean over multiple channels, or when input contains multiple subjects
%                        or trials)
%   cfg.maskstyle        = style used to masking, 'opacity', 'saturation', or 'outline' (default = 'opacity')
%                        'outline' can only be used with a logical cfg.maskparameter
%                        use 'saturation' or 'outline' when saving to vector-format (like *.eps) to avoid all sorts of image-problems
%   cfg.maskalpha        = alpha value between 0 (transparent) and 1 (opaque) used for masking areas dictated by cfg.maskparameter (default = 1)
%                        (will be ignored in case of numeric cfg.maskparameter or if cfg.maskstyle = 'outline')
%   cfg.masknans         = 'yes' or 'no' (default = 'yes')
%   cfg.xlim             = 'maxmin', 'maxabs', 'zeromax', 'minzero', or [xmin xmax] (default = 'maxmin')
%   cfg.ylim             = 'maxmin', 'maxabs', 'zeromax', 'minzero', or [ymin ymax] (default = 'maxmin')
%   cfg.zlim             = plotting limits for color dimension, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.gradscale        = number, scaling to apply to the MEG gradiometer channels prior to display
%   cfg.magscale         = number, scaling to apply to the MEG magnetometer channels prior to display
%   cfg.channel          = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.refchannel       = name of reference channel for visualising connectivity, can be 'gui'
%   cfg.baseline         = 'yes', 'no' or [time1 time2] (default = 'no'), see FT_FREQBASELINE
%   cfg.baselinetype     = 'absolute', 'relative', 'relchange', 'normchange', 'db' or 'zscore' (default = 'absolute')
%   cfg.trials           = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.box              = 'yes', 'no' (default = 'no' if maskparameter given default = 'yes')
%                          Draw a box around each graph
%   cfg.hotkeys          = enables hotkeys (up/down arrows) for dynamic colorbar adjustment
%   cfg.colorbar         = 'yes', 'no' (default = 'no')
%   cfg.colorbartext     = string indicating the text next to colorbar
%   cfg.colormap         = string, or Nx3 matrix, see FT_COLORMAP
%   cfg.showlabels       = 'yes', 'no' (default = 'no')
%   cfg.showoutline      = 'yes', 'no' (default = 'no')
%   cfg.showscale        = 'yes', 'no' (default = 'yes')
%   cfg.showcomment      = 'yes', 'no' (default = 'yes')
%   cfg.comment          = string of text (default = date + limits)
%                          Add 'comment' to graph (according to COMNT in the layout)
%   cfg.limittext        = add user-defined text instead of cfg.comment, (default = cfg.comment)
%   cfg.fontsize         = font size of comment and labels (if present) (default = 8)
%   cfg.fontweight       = font weight of comment and labels (if present)
%   cfg.interactive      = Interactive plot 'yes' or 'no' (default = 'yes')
%                          In a interactive plot you can select areas and produce a new
%                          interactive plot when a selected area is clicked. Multiple areas
%                          can be selected by holding down the SHIFT key.
%   cfg.renderer         = 'painters', 'zbuffer', ' opengl' or 'none' (default = [])
%   cfg.directionality   = '', 'inflow' or 'outflow' specifies for
%                          connectivity measures whether the inflow into a
%                          node, or the outflow from a node is plotted. The
%                          (default) behavior of this option depends on the dimor
%                          of the input data (see below).
%   cfg.layout           = specify the channel layout for plotting using one of
%                          the supported ways (see below).
%
% The following options for the scaling of the EEG, EOG, ECG, EMG, MEG and NIRS channels
% is optional and can be used to bring the absolute numbers of the different
% channel types in the same range (e.g. fT and uV). The channel types are determined
% from the input data using FT_CHANNELSELECTION.
%   cfg.eegscale         = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale         = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale         = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale         = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale         = number, scaling to apply to the MEG channels prior to display
%   cfg.gradscale        = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale         = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.nirsscale        = number, scaling to apply to the NIRS channels prior to display
%   cfg.mychanscale      = number, scaling to apply to the channels specified in cfg.mychan
%   cfg.mychan           = Nx1 cell-array with selection of channels
%   cfg.chanscale        = Nx1 vector with scaling factors, one per channel specified in cfg.channel
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
% The layout defines how the channels are arranged and what the size of each
% subplot is. You can specify the layout in a variety of ways:
%  - you can provide a pre-computed layout structure, see FT_PREPARE_LAYOUT
%  - you can give the name of an ASCII layout file with extension *.lay
%  - you can give the name of an electrode file
%  - you can give an electrode definition, i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% If you do not specify any of these and the data structure contains an
% electrode or gradiometer structure (common for MEG data, since the header
% of the MEG datafile contains the gradiometer information), that will be
% used for creating a layout. If you want to have more fine-grained control
% over the layout of the subplots, you should create your own layout file.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure. For this particular function, the
% data should be provided as a cell-array.
%
% See also:
%   FT_MULTIPLOTER, FT_SINGLEPLOTER, FT_SINGLEPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR,
%   FT_PREPARE_LAYOUT

% Undocumented local options:
% cfg.channel
% cfg.layoutname
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
ft_preamble loadvar data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'freq');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels', 'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamed',    {'cohrefchannel', 'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'matrixside', 'directionality'});
cfg = ft_checkconfig(cfg, 'renamed',    {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedback', 'inflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim', 'absmax', 'maxabs'});
cfg = ft_checkconfig(cfg, 'unused',     {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'newfigure', 'figure'});

% set the defaults
cfg.parameter      = ft_getopt(cfg, 'parameter',    'powspctrm');
cfg.baseline       = ft_getopt(cfg, 'baseline',     'no');
cfg.baselinetype   = ft_getopt(cfg, 'baselinetype', 'absolute');
cfg.trials         = ft_getopt(cfg, 'trials', 'all', 1);
cfg.xlim           = ft_getopt(cfg, 'xlim',         'maxmin');
cfg.ylim           = ft_getopt(cfg, 'ylim',         'maxmin');
cfg.zlim           = ft_getopt(cfg, 'zlim',         'maxmin');
cfg.colorbar       = ft_getopt(cfg, 'colorbar',     'no');
cfg.colormap       = ft_getopt(cfg, 'colormap',      'default');
cfg.colorbartext   = ft_getopt(cfg, 'colorbartext', '');
cfg.comment        = ft_getopt(cfg, 'comment',       date);
cfg.limittext      = ft_getopt(cfg, 'limittext',    'default');
cfg.showlabels     = ft_getopt(cfg, 'showlabels',   'no');
cfg.showoutline    = ft_getopt(cfg, 'showoutline',  'no');
cfg.showscale      = ft_getopt(cfg, 'showscale',    'yes');
cfg.showcomment    = ft_getopt(cfg, 'showcomment',  'yes');
cfg.channel        = ft_getopt(cfg, 'channel',      'all');
cfg.fontsize       = ft_getopt(cfg, 'fontsize',      8);
cfg.fontweight     = ft_getopt(cfg, 'fontweight');
cfg.interactive    = ft_getopt(cfg, 'interactive',  'yes');
cfg.hotkeys        = ft_getopt(cfg, 'hotkeys',      'yes');
cfg.orient         = ft_getopt(cfg, 'orient',       'landscape');
cfg.maskalpha      = ft_getopt(cfg, 'maskalpha',     1);
cfg.masknans       = ft_getopt(cfg, 'masknans',     'yes');
cfg.maskparameter  = ft_getopt(cfg, 'maskparameter');
cfg.maskstyle      = ft_getopt(cfg, 'maskstyle',    'opacity');
cfg.directionality = ft_getopt(cfg, 'directionality', '');
cfg.figurename     = ft_getopt(cfg, 'figurename');
cfg.commentpos     = ft_getopt(cfg, 'commentpos',   'layout');
cfg.scalepos       = ft_getopt(cfg, 'scalepos',     'layout');
cfg.renderer       = ft_getopt(cfg, 'renderer'); % let MATLAB decide on the default

if ~isfield(cfg, 'box')
  if ~isempty(cfg.maskparameter)
    cfg.box = 'yes';
  else
    cfg.box = 'no';
  end
end

% check if the colormap is in the proper format
if ~isequal(cfg.colormap, 'default')
  if ischar(cfg.colormap)
    cfg.colormap = ft_colormap(cfg.colormap);
  elseif iscell(cfg.colormap)
    cfg.colormap = ft_colormap(cfg.colormap{:});
  elseif isnumeric(cfg.colormap) && size(cfg.colormap,2)~=3
    ft_error('cfg.colormap must be Nx3');
  end
  % the actual colormap will be set below
end

% this is needed for the figure title
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

hastime = isfield(data, 'time');
hasfreq = isfield(data, 'freq');

assert((hastime && hasfreq), 'please use ft_multiplotER for time-only or frequency-only data');

xparam = ft_getopt(cfg, 'xparam', 'time');
yparam = ft_getopt(cfg, 'yparam', 'freq');

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
  tmpcfg = removefields(cfg, {'inputfile', 'reproducescript'});
  % keep mask-parameter if it is set
  if ~isempty(cfg.maskparameter)
    tempmask = data.(cfg.maskparameter);
  end
  data = ft_freqbaseline(tmpcfg, data);
  % put mask-parameter back if it is set
  if ~isempty(cfg.maskparameter)
    data.(cfg.maskparameter) = tempmask;
  end
end

% channels SHOULD be selected here, as no interactive action produces a new multiplot
tmpcfg = keepfields(cfg, {'channel', 'trials', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
if hasrpt
  tmpcfg.avgoverrpt = 'yes';
else
  tmpcfg.avgoverrpt = 'no';
end
tmpvar = data;
[data] = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

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
% refchannel can relate to connectivity or to an EEG reference.

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
  % convert the bivariate data to univariate and call this plotting function with univariate input
  cfg.originalfunction = 'ft_multiplotTFR';
  cfg.trials = 'all'; % trial selection has been taken care off
  bivariate_common(cfg, data);
  return
end

% Apply channel-type specific scaling
fn = fieldnames(cfg);
fn = setdiff(fn, {'skipscale', 'showscale', 'gridscale'}); % these are for the layout and plotting, not for CHANSCALE_COMMON
fn = fn(endsWith(fn, 'scale') | startsWith(fn, 'mychan') | strcmp(fn, 'channel') | strcmp(fn, 'parameter'));
tmpcfg = keepfields(cfg, fn);
if ~isempty(tmpcfg)
  data = chanscale_common(tmpcfg, data);
  % remove the scaling fields from the configuration, to prevent them from being called again in interactive mode
  % but keep the parameter and channel field
  cfg = removefields(cfg, setdiff(fn, {'parameter', 'channel'}));
else
  % do nothing
end


%% Section 3: select the data to be plotted and determine min/max range

% Read or create the layout that will be used for plotting
tmpcfg = keepfields(cfg, {'layout', 'channel', 'rows', 'columns', 'commentpos', 'skipcomnt', 'scalepos', 'skipscale', 'projection', 'viewpoint', 'rotate', 'width', 'height', 'elec', 'grad', 'opto', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
cfg.layout = ft_prepare_layout(tmpcfg, data);

% Take the subselection of channels that is contained in the layout, this is the same in all datasets
[selchan, sellay] = match_str(data.label, cfg.layout.label);

% Get physical min/max range of x, i.e. time
if ~isnumeric(cfg.xlim)
  xmin = nanmin(data.(xparam));
  xmax = nanmax(data.(xparam));
  switch cfg.xlim
    case 'maxmin'
      % keep them as they are
    case 'maxabs'
      xmax = max(abs(xmax), abs(xmin));
      xmin = -xmax;
    case 'zeromax'
      xmin = 0;
    case 'minzero'
      xmax = 0;
    otherwise
      ft_error('invalid specification of cfg.xlim');
  end % switch
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Get the index of the nearest bin, this is the same in all datasets
xminindx = nearest(data.(xparam), xmin);
xmaxindx = nearest(data.(xparam), xmax);
xmin = data.(xparam)(xminindx);
xmax = data.(xparam)(xmaxindx);
selx = xminindx:xmaxindx;
xval = data.(xparam)(selx);

% Get physical min/max range of y, i.e. freq
if ~isnumeric(cfg.ylim)
  ymin = nanmin(data.(yparam));
  ymax = nanmax(data.(yparam));
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
  end % switch
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Get the index of the nearest bin, this is the same in all datasets
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
[c, ia, ib] = intersect({'chan', yparam, xparam}, dimtok, 'stable');
datamatrix = permute(datamatrix, ib);
datamatrix = datamatrix(selchan, sely, selx);

if ~isempty(cfg.maskparameter)
  % one value for each channel-freq-time point
  maskmatrix = data.(cfg.maskparameter)(selchan, sely, selx);
  if islogical(maskmatrix) && any(strcmp(cfg.maskstyle, {'saturation', 'opacity'}))
    maskmatrix = double(maskmatrix);
    maskmatrix(~maskmatrix) = cfg.maskalpha;
  elseif isnumeric(maskmatrix)
    if strcmp(cfg.maskstyle, 'outline')
      error('Outline masking with a numeric cfg.maskparameter is not supported. Please use a logical mask instead.')
    end
    if cfg.maskalpha ~= 1
      ft_warning('Using field "%s" for masking, cfg.maskalpha is ignored.', cfg.maskparameter)
    end
    % scale mask between 0 and 1
    minval = min(maskmatrix(:));
    maxval = max(maskmatrix(:));
    maskmatrix = (maskmatrix - minval) / (maxval-minval);
  end
else
  % create an Nx0x0 matrix
  maskmatrix = zeros(length(selchan), 0, 0);
end

chanX      = cfg.layout.pos(sellay, 1);
chanY      = cfg.layout.pos(sellay, 2);
chanWidth  = cfg.layout.width(sellay);
chanHeight = cfg.layout.height(sellay);

%% Section 4: do the actual plotting

% open a new figure, or add it to the existing one
% note that in general adding a TFR to an existing one does not make sense, since they will overlap
open_figure(keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'}));
hold on

% Get physical z-axis range (color axis):
if strcmp(cfg.zlim, 'maxmin')
  zmin = nanmin(datamatrix(:));
  zmax = nanmax(datamatrix(:));
elseif strcmp(cfg.zlim, 'maxabs')
  zmin = -nanmax(abs(datamatrix(:)));
  zmax =  nanmax(abs(datamatrix(:)));
elseif strcmp(cfg.zlim, 'zeromax')
  zmin = 0;
  zmax = nanmax(datamatrix(:));
elseif strcmp(cfg.zlim, 'minzero')
  zmin = nanmin(datamatrix(:));
  zmax = 0;
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% Plot channels
for k=1:length(selchan)
  cdata = shiftdim(datamatrix(k, :, :));
  if ~isempty(cfg.maskparameter)
    mdata = shiftdim(maskmatrix(k, :, :));
  end
  
  % Draw plot (and mask Nan's with maskfield if requested)
  if isequal(cfg.masknans, 'yes') && isempty(cfg.maskparameter)
    nans_mask = ~isnan(cdata);
    mask = double(nans_mask);
    ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask, 'hpos', chanX(k), 'vpos', chanY(k), 'width', chanWidth(k), 'height', chanHeight(k))
  elseif isequal(cfg.masknans, 'yes') && ~isempty(cfg.maskparameter)
    nans_mask = ~isnan(cdata);
    mask = nans_mask .* mdata;
    mask = double(mask);
    ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask, 'hpos', chanX(k), 'vpos', chanY(k), 'width', chanWidth(k), 'height', chanHeight(k))
  elseif isequal(cfg.masknans, 'no') && ~isempty(cfg.maskparameter)
    mask = mdata;
    mask = double(mask);
    ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask, 'hpos', chanX(k), 'vpos', chanY(k), 'width', chanWidth(k), 'height', chanHeight(k))
  else
    ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'hpos', chanX(k), 'vpos', chanY(k), 'width', chanWidth(k), 'height', chanHeight(k))
  end
  
  % Currently the handle isn't being used below, this is here for possible use in the future
  h = findobj('tag', 'cip');
end % plot channels

% plot the layout, labels and outline
ft_plot_layout(cfg.layout, 'box', cfg.box, 'label', cfg.showlabels, 'outline', cfg.showoutline, 'point', 'no', 'mask', 'no', 'fontsize', cfg.fontsize, 'labelyoffset', 1.4*median(cfg.layout.height/2), 'labelalignh', 'center', 'chanindx', find(~ismember(cfg.layout.label, {'COMNT', 'SCALE'})) );

% apply the colormap
if ~isempty(cfg.colormap)
  set(gcf,  'colormap', cfg.colormap);
end

% Construct comment
if ~strcmp(cfg.comment, 'no')
  if ~strcmp(cfg.limittext, 'default')
    comment = cfg.limittext;
  else
    comment = cfg.comment;
    comment = sprintf('%0s\nxlim=[%.3g %.3g]', comment, xmin, xmax);
    comment = sprintf('%0s\nylim=[%.3g %.3g]', comment, ymin, ymax);
    comment = sprintf('%0s\nzlim=[%.3g %.3g]', comment, zmin, zmax);
  end
end

% Write comment
if strcmp(cfg.comment, 'no')
  comment_handle = [];
elseif strcmp(cfg.commentpos, 'title')
  comment_handle = title(comment, 'FontSize', cfg.fontsize);
elseif ~isempty(strcmp(cfg.layout.label, 'COMNT'))
  x_comment = cfg.layout.pos(strcmp(cfg.layout.label, 'COMNT'), 1);
  y_comment = cfg.layout.pos(strcmp(cfg.layout.label, 'COMNT'), 2);
  % 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',
  comment_handle = ft_plot_text(x_comment, y_comment, comment, 'FontSize', cfg.fontsize, 'FontWeight', cfg.fontweight);
else
  comment_handle = [];
end

% show scale
if istrue(cfg.showscale)
  k = cellstrmatch('SCALE', cfg.layout.label);
  if ~isempty(k)
    % Get average cdata across channels:
    cdata = shiftdim(mean(datamatrix, 1));
    
    % Draw plot (and mask Nan's with maskfield if requested)
    if isequal(cfg.masknans, 'yes') && isempty(cfg.maskparameter)
      mask = ~isnan(cdata);
      mask = double(mask);
      ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask, 'hpos', cfg.layout.pos(k, 1), 'vpos', cfg.layout.pos(k, 2), 'width', cfg.layout.width(k), 'height', cfg.layout.height(k))
    elseif isequal(cfg.masknans, 'yes') && ~isempty(cfg.maskparameter)
      mask = ~isnan(cdata);
      mask = mask .* mdata;
      mask = double(mask);
      ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask, 'hpos', cfg.layout.pos(k, 1), 'vpos', cfg.layout.pos(k, 2), 'width', cfg.layout.width(k), 'height', cfg.layout.height(k))
    elseif isequal(cfg.masknans, 'no') && ~isempty(cfg.maskparameter)
      mask = mdata;
      mask = double(mask);
      ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask, 'hpos', cfg.layout.pos(k, 1), 'vpos', cfg.layout.pos(k, 2), 'width', cfg.layout.width(k), 'height', cfg.layout.height(k))
    else
      ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'hpos', cfg.layout.pos(k, 1), 'vpos', cfg.layout.pos(k, 2), 'width', cfg.layout.width(k), 'height', cfg.layout.height(k))
    end
  end
end % show scale

% show colorbar
if isfield(cfg, 'colorbar') && (strcmp(cfg.colorbar, 'yes'))
  c = colorbar;
  ylabel(c, cfg.colorbartext);
end

% Set colour axis
caxis([zmin zmax]);
if strcmp('yes', cfg.hotkeys)
  %  Attach data and cfg to figure and attach a key listener to the figure
  set(gcf, 'KeyPressFcn', {@key_sub, zmin, zmax})
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

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  % add the cfg/data/channel information to the figure under identifier linked to this axis
  ident                 = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
  set(gca,'tag',ident);
  info                  = guidata(gcf);
  info.(ident).x        = cfg.layout.pos(:, 1);
  info.(ident).y        = cfg.layout.pos(:, 2);
  info.(ident).label    = cfg.layout.label;
  info.(ident).dataname = dataname;
  info.(ident).cfg      = cfg;
  info.(ident).data     = data;
  info.(ident).commenth = comment_handle;
  guidata(gcf, info);
  set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR}, 'event', 'WindowButtonMotionFcn'});
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous data
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
function l = cellstrmatch(str, strlist)
l = [];
for k=1:length(strlist)
  if strcmp(char(str), char(strlist(k)))
    l = [l k];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotTFR(label, varargin)
% fetch cfg/data based on axis indentifier given as tag
ident = get(gca,'tag');
info  = guidata(gcf);
cfg   = info.(ident).cfg;
data  = info.(ident).data;
if ~isempty(label)
  cfg = removefields(cfg, 'inputfile');   % the reading has already been done and varargin contains the data
  cfg.baseline = 'no';                    % make sure the next function does not apply a baseline correction again
  cfg.channel = label;
  cfg.dataname = info.(ident).dataname;   % put data name in here, this cannot be resolved by other means
  cfg.trials = 'all';                     % trial selection has already been taken care of
  fprintf('selected cfg.channel = {%s}\n', join_str(', ', cfg.channel));
  % ensure that the new figure appears at the same position
  cfg.figure = 'yes';
  cfg.position = get(gcf, 'Position');
  ft_singleplotTFR(cfg, data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata, varargin)
ident       = get(gca, 'tag');
info        = guidata(gcf);

climits = caxis;
incr_c  = abs(climits(2) - climits(1)) /10;

newz = climits;
if length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:}, 'control')
  % TRANSLATE by 10%
  switch eventdata.Key
    case 'pageup'
      newz = [climits(1)+incr_c climits(2)+incr_c];
    case 'pagedown'
      newz = [climits(1)-incr_c climits(2)-incr_c];
  end % switch
else
  % ZOOM by 10%
  switch eventdata.Key
    case 'pageup'
      newz = [climits(1)-incr_c climits(2)+incr_c];
    case 'pagedown'
      newz = [climits(1)+incr_c climits(2)-incr_c];
    case 'm'
      newz = [varargin{1} varargin{2}];
  end % switch
end % if

% update the color axis
caxis(newz);

if ~isempty(ident) && isfield(info.(ident), 'commenth') && ~isempty(info.(ident).commenth)
  commentstr = get(info.(ident).commenth, 'string');
  sel        = contains(commentstr, 'zlim');
  if any(sel)
    commentstr{sel} = sprintf('%0s=[%.3g %.3g]', 'zlim', newz(1), newz(2));
    set(info.(ident).commenth, 'string', commentstr);
  end
end
