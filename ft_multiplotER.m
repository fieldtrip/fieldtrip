function [cfg] = ft_multiplotER(cfg, varargin)

% FT_MULTIPLOTER plots the event-related potentials, event-related fields
% or oscillatory activity (power or coherence) versus frequency. Multiple
% datasets can be overlayed.  The plots are arranged according to their
% location specified in the layout.
%
% Use as
%   ft_multiplotER(cfg, data)
% or
%   ft_multiplotER(cfg, data, data2, ..., dataN)
%
% The data can be an ERP/ERF produced by FT_TIMELOCKANALYSIS, a powerspectrum
% produced by FT_FREQANALYSIS or a coherencespectrum produced by FT_FREQDESCRIPTIVES.
% If you specify multiple datasets they must contain the same channels, etc.
%
% The configuration can have the following parameters:
%   cfg.parameter     = field to be plotted on y-axis (default depends on data.dimord)
%                       'avg', 'powspctrm' or 'cohspctrm'
%   cfg.maskparameter = field in the first dataset to be used for marking significant data
%   cfg.maskstyle     = style used for masking of data, 'box', 'thickness' or 'saturation' (default = 'box')
%   cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim          = 'maxmin', 'maxabs', 'zeromax', 'minzero', or [ymin ymax] (default = 'maxmin')
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.refchannel    = name of reference channel for visualising connectivity, can be 'gui'
%   cfg.baseline      = 'yes', 'no' or [time1 time2] (default = 'no'), see FT_TIMELOCKBASELINE or FT_FREQBASELINE
%   cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.axes          = 'yes', 'no' (default = 'yes')
%                       Draw x- and y-axes for each graph
%   cfg.box           = 'yes', 'no' (default = 'no')
%                       Draw a box around each graph
%   cfg.comment       = string of text (default = date + colors)
%                       Add 'comment' to graph (according to COMNT in the layout)
%   cfg.showlabels    = 'yes', 'no' (default = 'no')
%   cfg.showoutline   = 'yes', 'no' (default = 'no')
%   cfg.fontsize      = font size of comment and labels (if present) (default = 8)
%   cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'yes')
%                       In a interactive plot you can select areas and produce a new
%                       interactive plot when a selected area is clicked. Multiple areas
%                       can be selected by holding down the SHIFT key.
%   cfg.renderer      = 'painters', 'zbuffer', ' opengl' or 'none' (default = [])
%   cfg.linestyle     = linestyle/marker type, see options of the PLOT function (default = '-')
%                       can be a single style for all datasets, or a cell-array containing one style for each dataset
%   cfg.linewidth     = linewidth in points (default = 0.5)
%   cfg.graphcolor    = color(s) used for plotting the dataset(s) (default = 'brgkywrgbkywrgbkywrgbkyw')
%                       alternatively, colors can be specified as Nx3 matrix of RGB values
%   cfg.directionality = '', 'inflow' or 'outflow' specifies for
%                       connectivity measures whether the inflow into a
%                       node, or the outflow from a node is plotted. The
%                       (default) behavior of this option depends on the dimor
%                       of the input data (see below).
%   cfg.layout        = specify the channel layout for plotting using one of
%                       the supported ways (see below).
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
% data should be provided as a cell array.
%
% See also FT_MULTIPLOTTFR, FT_SINGLEPLOTER, FT_SINGLEPLOTTFR, FT_TOPOPLOTER,
% FT_TOPOPLOTTFR, FT_PREPARE_LAYOUT

% Undocumented local options:
% cfg.layoutname
% cfg.preproc
% cfg.orient = landscape/portrait

% Copyright (C) 2003-2006, Ole Jensen
% Copyright (C) 2007-2011, Roemer van der Meij & Jan-Mathijs Schoffelen
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

for i=1:length(varargin)
  % check if the input data is valid for this function
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'timelock', 'freq'});
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim', 'absmax', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedback', 'inflow'});
cfg = ft_checkconfig(cfg, 'renamed', {'matrixside',  'directionality'});
cfg = ft_checkconfig(cfg, 'renamed', {'cohrefchannel', 'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed', {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'renamed', {'hlim', 'xlim'});
cfg = ft_checkconfig(cfg, 'renamed', {'vlim', 'ylim'});
cfg = ft_checkconfig(cfg, 'deprecated', {'xparam'});
cfg = ft_checkconfig(cfg, 'unused',  {'cohtargetchannel'});

% set the defaults
cfg.baseline       = ft_getopt(cfg, 'baseline', 'no');
cfg.trials         = ft_getopt(cfg, 'trials', 'all', 1);
cfg.xlim           = ft_getopt(cfg, 'xlim', 'maxmin');
cfg.ylim           = ft_getopt(cfg, 'ylim', 'maxmin');
cfg.comment        = ft_getopt(cfg, 'comment', strcat([date '\n']));
cfg.axes           = ft_getopt(cfg, 'axes', 'yes');
cfg.showlabels     = ft_getopt(cfg, 'showlabels', 'no');
cfg.showoutline    = ft_getopt(cfg, 'showoutline', 'no');
cfg.box            = ft_getopt(cfg, 'box', 'no');
cfg.fontsize       = ft_getopt(cfg, 'fontsize', 8);
cfg.graphcolor     = ft_getopt(cfg, 'graphcolor', 'brgkywrgbkywrgbkywrgbkyw');
cfg.interactive    = ft_getopt(cfg, 'interactive', 'yes');
cfg.renderer       = ft_getopt(cfg, 'renderer'); % let MATLAB decide on default
cfg.orient         = ft_getopt(cfg, 'orient', 'landscape');
cfg.maskparameter  = ft_getopt(cfg, 'maskparameter');
cfg.linestyle      = ft_getopt(cfg, 'linestyle', '-');
cfg.linewidth      = ft_getopt(cfg, 'linewidth', 0.5);
cfg.maskstyle      = ft_getopt(cfg, 'maskstyle', 'box');
cfg.channel        = ft_getopt(cfg, 'channel', 'all');
cfg.directionality = ft_getopt(cfg, 'directionality', '');
cfg.figurename     = ft_getopt(cfg, 'figurename');
cfg.preproc        = ft_getopt(cfg, 'preproc');
cfg.tolerance      = ft_getopt(cfg, 'tolerance', 1e-5);
cfg.frequency      = ft_getopt(cfg, 'frequency', 'all'); % needed for frequency selection with TFR data
cfg.latency        = ft_getopt(cfg, 'latency', 'all'); % needed for latency selection with TFR data, FIXME, probably not used

if numel(findobj(gcf, 'type', 'axes', '-not', 'tag', 'ft-colorbar')) > 1 && strcmp(cfg.interactive, 'yes')
  warning('using cfg.interactive = ''yes'' in subplots is not supported, setting cfg.interactive = ''no''')
  cfg.interactive = 'no';
end

Ndata = length(varargin);

for i=1:Ndata
  dtype{i}   = ft_datatype(varargin{i});
  hastime(i) = ~isempty(strfind(varargin{i}.dimord, 'time'));
  hasfreq(i) = ~isempty(strfind(varargin{i}.dimord, 'freq'));
end

% check if the input has consistent datatypes
if ~all(strcmp(dtype, dtype{1})) || ~all(hastime==hastime(1)) || ~all(hasfreq==hasfreq(1))
  error('different datatypes are not allowed as input');
end
dtype   = dtype{1};
hastime = hastime(1);
hasfreq = hasfreq(1);

% ensure that all inputs are sufficiently consistent
if hastime && ~checktime(varargin{:}, 'identical', cfg.tolerance);
  error('this function requires identical time axes for all input structures');
end
if hasfreq && ~checkfreq(varargin{:}, 'identical', cfg.tolerance);
  error('this function requires identical frequency axes for all input structures');
end

%FIXME rename directionality and refchannel in more meaningful options
if ischar(cfg.graphcolor)
  GRAPHCOLOR = ['k' cfg.graphcolor];
elseif isnumeric(cfg.graphcolor)
  GRAPHCOLOR = [0 0 0; cfg.graphcolor];
end

% check for linestyle being a cell-array, check it's length, and lengthen it if does not have enough styles in it
if ischar(cfg.linestyle)
  cfg.linestyle = {cfg.linestyle};
end

if Ndata>1
  if (length(cfg.linestyle) < Ndata ) && (length(cfg.linestyle) > 1)
    error('either specify cfg.linestyle as a cell-array with one cell for each dataset, or only specify one linestyle')
  elseif (length(cfg.linestyle) < Ndata ) && (length(cfg.linestyle) == 1)
    tmpstyle = cfg.linestyle{1};
    cfg.linestyle = cell(Ndata , 1);
    for idataset = 1:Ndata
      cfg.linestyle{idataset} = tmpstyle;
    end
  end
end

% % interactive plotting is not allowed with more than 1 input
% if numel(varargin)>1 && strcmp(cfg.interactive, 'yes')
%   error('interactive plotting is not supported with more than 1 input data set');
% end

dimord = varargin{1}.dimord;
dimtok = tokenize(dimord, '_');

% ensure that the preproc specific options are located in the cfg.preproc
% substructure, but also ensure that the field 'refchannel' is present at the
% highest level in the structure. This is a little hack by JM because the field
% refchannel can also refer to the plotting of a connectivity metric. Also,
% the freq2raw conversion does not work at all in the call to ft_preprocessing.
% Therefore, for now, the preprocessing will not be done when there is freq
% data in the input. A more generic solution should be considered.

if isfield(cfg, 'refchannel'), refchannelincfg = cfg.refchannel; end
if ~any(strcmp({'freq', 'freqmvar'}, dtype)),
  cfg = ft_checkconfig(cfg, 'createsubcfg', {'preproc'});
end
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

for i=1:Ndata
  % this is needed for correct treatment of GRAPHCOLOR later on
  if nargin>1,
    if ~isempty(inputname(i+1))
      iname{i+1} = inputname(i+1);
    else
      iname{i+1} = ['input', num2str(i, '%02d')];
    end
  else
    iname{i+1} = cfg.inputfile{i};
  end
end

% Set x/y/parameter defaults according to datatype and dimord
switch dtype
  case 'timelock'
    xparam = 'time';
    yparam = '';
    cfg.parameter = ft_getopt(cfg, 'parameter', 'avg');
  case 'freq'
    if any(ismember(dimtok, 'time'))
      xparam = 'time';
      yparam = 'freq';
      cfg.parameter = ft_getopt(cfg, 'parameter', 'powspctrm');
    else
      xparam = 'freq';
      yparam = '';
      cfg.parameter = ft_getopt(cfg, 'parameter', 'powspctrm');
    end
  case 'comp'
    % not supported
  otherwise
    % not supported
end

% user specified own fields, but no yparam (which is not asked in help)
if exist('xparam', 'var') && isfield(cfg, 'parameter') && ~exist('yparam', 'var')
  yparam = '';
end

if isfield(cfg, 'channel') && isfield(varargin{1}, 'label')
  cfg.channel = ft_channelselection(cfg.channel, varargin{1}.label);
elseif isfield(cfg, 'channel') && isfield(varargin{1}, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(varargin{1}.labelcmb(:)));
end

% perform channel selection, unless in the other plotting functions this
% can always be done because ft_multiplotER is the entry point into the
% interactive stream, but will not be revisited
if isfield(varargin{1}, 'label')
  % only do the channel selection when it can actually be done,
  % i.e. when the data are bivariate ft_selectdata will crash, moreover
  % the bivariate case is handled below
  tmpcfg = keepfields(cfg, 'channel');
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

  clear tmpvar tmpcfg
end

if isfield(varargin{1}, 'label') % && strcmp(cfg.interactive, 'no')
  selchannel = ft_channelselection(cfg.channel, varargin{1}.label);
elseif isfield(varargin{1}, 'labelcmb') % && strcmp(cfg.interactive, 'no')
  selchannel = ft_channelselection(cfg.channel, unique(varargin{1}.labelcmb(:)));
end

% check whether rpt/subj is present and remove if necessary
% FIXME this should be implemented with avgoverpt in ft_selectdata
hasrpt = sum(ismember(dimtok, {'rpt' 'subj'}));
if strcmp(dtype, 'timelock') && hasrpt,
  tmpcfg = [];

  % disable hashing of input data (speeds up things)
  tmpcfg.trackcallinfo = 'no';

  tmpcfg.trials = cfg.trials;
  for i=1:Ndata
    % save mask (timelockanalysis will remove it)
    if ~isempty(cfg.maskparameter)
      tmpmask = varargin{i}.(cfg.maskparameter);
    end
    varargin{i} = ft_timelockanalysis(tmpcfg, varargin{i});
    if ~strcmp(cfg.parameter, 'avg')
      % rename avg back into its original parameter name
      varargin{i}.(cfg.parameter) = varargin{i}.avg;
      varargin{i} = rmfield(varargin{i}, 'avg');
    end

    % put back mask
    if ~isempty(cfg.maskparameter)
      varargin{i}.(cfg.maskparameter) = tmpmask;
    end
  end
  dimord        = varargin{1}.dimord;
  dimtok        = tokenize(dimord, '_');

elseif strcmp(dtype, 'freq') && hasrpt,
  % this also deals with fourier-spectra in the input
  % or with multiple subjects in a frequency domain stat-structure
  % on the fly computation of coherence spectrum is not supported
  for i=1:Ndata
    if isfield(varargin{i}, 'crsspctrm'),
      varargin{i} = rmfield(varargin{i}, 'crsspctrm');
    end
  end

  tmpcfg           = [];
  tmpcfg.trials    = cfg.trials;
  tmpcfg.jackknife = 'no';
  for i=1:Ndata
    if isfield(cfg, 'parameter') && ~strcmp(cfg.parameter, 'powspctrm')
      % freqdesctiptives will only work on the powspctrm field
      % hence a temporary copy of the data is needed
      tempdata.dimord    = varargin{i}.dimord;
      tempdata.freq      = varargin{i}.freq;
      tempdata.label     = varargin{i}.label;
      tempdata.powspctrm = varargin{i}.(cfg.parameter);
      if isfield(varargin{i}, 'cfg') tempdata.cfg = varargin{i}.cfg; end
      tempdata           = ft_freqdescriptives(tmpcfg, tempdata);
      varargin{i}.(cfg.parameter)  = tempdata.powspctrm;
      clear tempdata
    else
      varargin{i} = ft_freqdescriptives(tmpcfg, varargin{i});
    end
  end
  dimord = varargin{1}.dimord;
  dimtok = tokenize(dimord, '_');
end

% % Read or create the layout that will be used for plotting
cla
lay = ft_prepare_layout(cfg, varargin{1});
cfg.layout = lay;

% plot layout
boxflg     = istrue(cfg.box);
labelflg   = false; % channel labels are plotted further down using ft_plot_vector
outlineflg = istrue(cfg.showoutline);
ft_plot_lay(lay, 'box', boxflg, 'label', labelflg, 'outline', outlineflg, 'point', 'no', 'mask', 'no');

% Apply baseline correction
if ~strcmp(cfg.baseline, 'no')
  for i=1:Ndata
    if strcmp(dtype, 'timelock') && strcmp(xparam, 'time')
      varargin{i} = ft_timelockbaseline(cfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(xparam, 'time')
      varargin{i} = ft_freqbaseline(cfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(xparam, 'freq')
      error('Baseline correction is not supported for spectra without a time dimension');
    else
      warning('Baseline correction not applied, please set xparam');
    end
  end
end

% Handle the bivariate case

% Check for bivariate metric with 'chan_chan' in the dimord
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;

% Check for bivariate metric with a labelcmb
haslabelcmb = isfield(varargin{1}, 'labelcmb');

if (isfull || haslabelcmb) && (isfield(varargin{1}, cfg.parameter) && ~strcmp(cfg.parameter, 'powspctrm'))
  % A reference channel is required:
  if ~isfield(cfg, 'refchannel')
    error('no reference channel is specified');
  end

  % check for refchannel being part of selection
  if ~strcmp(cfg.refchannel, 'gui')
    if haslabelcmb
      cfg.refchannel = ft_channelselection(cfg.refchannel, unique(varargin{1}.labelcmb(:)));
    else
      cfg.refchannel = ft_channelselection(cfg.refchannel, varargin{1}.label);
    end
    if (isfull      && ~any(ismember(varargin{1}.label, cfg.refchannel))) || ...
        (haslabelcmb && ~any(ismember(varargin{1}.labelcmb(:), cfg.refchannel)))
      error('cfg.refchannel is a not present in the (selected) channels)')
    end
  end

  % Interactively select the reference channel
  if strcmp(cfg.refchannel, 'gui')
    % Open a single figure with the channel layout, the user can click on a reference channel
    h = clf;
    ft_plot_lay(lay, 'box', false);
    title('Select the reference channel by dragging a selection window, more than 1 channel can be selected...');
    % add the channel information to the figure
    info       = guidata(gcf);
    info.x     = lay.pos(:, 1);
    info.y     = lay.pos(:, 2);
    info.label = lay.label;
    guidata(h, info);
    %set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'callback', {@select_topoplotER, cfg, data}});
    set(gcf, 'WindowButtonUpFcn',  {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotER, cfg, varargin{1}}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotER, cfg, varargin{1}}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotER, cfg, varargin{1}}, 'event', 'WindowButtonMotionFcn'});
    return
  end

  for i=1:Ndata
    if ~isfull,
      % Convert 2-dimensional channel matrix to a single dimension:
      if isempty(cfg.directionality)
        sel1 = find(strcmp(cfg.refchannel, varargin{i}.labelcmb(:, 2)));
        sel2 = find(strcmp(cfg.refchannel, varargin{i}.labelcmb(:, 1)));
      elseif strcmp(cfg.directionality, 'outflow')
        sel1 = [];
        sel2 = find(strcmp(cfg.refchannel, varargin{i}.labelcmb(:, 1)));
      elseif strcmp(cfg.directionality, 'inflow')
        sel1 = find(strcmp(cfg.refchannel, varargin{i}.labelcmb(:, 2)));
        sel2 = [];
      end
      fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.parameter);
      if length(sel1)+length(sel2)==0
        error('there are no channels selected for plotting: you may need to look at the specification of cfg.directionality');
      end
      varargin{i}.(cfg.parameter) = varargin{i}.(cfg.parameter)([sel1;sel2], :, :);
      varargin{i}.label     = [varargin{i}.labelcmb(sel1, 1);varargin{i}.labelcmb(sel2, 2)];
      varargin{i}.labelcmb  = varargin{i}.labelcmb([sel1;sel2], :);
      %varargin{i}           = rmfield(varargin{i}, 'labelcmb');
    else
      % General case
      sel               = match_str(varargin{i}.label, cfg.refchannel);
      siz               = [size(varargin{i}.(cfg.parameter)) 1];
      if strcmp(cfg.directionality, 'inflow') || isempty(cfg.directionality)
        %the interpretation of 'inflow' and 'outflow' depend on
        %the definition in the bivariate representation of the data
        %in FieldTrip the row index 'causes' the column index channel
        %data.(cfg.parameter) = reshape(mean(data.(cfg.parameter)(:, sel, :), 2), [siz(1) 1 siz(3:end)]);
        sel1 = 1:siz(1);
        sel2 = sel;
        meandir = 2;
      elseif strcmp(cfg.directionality, 'outflow')
        %data.(cfg.parameter) = reshape(mean(data.(cfg.parameter)(sel, :, :), 1), [siz(1) 1 siz(3:end)]);
        sel1 = sel;
        sel2 = 1:siz(1);
        meandir = 1;

      elseif strcmp(cfg.directionality, 'ff-fd')
        error('cfg.directionality = ''ff-fd'' is not supported anymore, you have to manually subtract the two before the call to ft_multiplotER');
      elseif strcmp(cfg.directionality, 'fd-ff')
        error('cfg.directionality = ''fd-ff'' is not supported anymore, you have to manually subtract the two before the call to ft_multiplotER');
      end %if directionality
    end %if ~isfull
  end %for i
end %handle the bivariate data

% Get physical min/max range of x
if strcmp(cfg.xlim, 'maxmin')
  % Find maxmin throughout all varargins:
  xmin = [];
  xmax = [];
  for i=1:length(varargin)
    xmin = min([xmin varargin{i}.(xparam)]);
    xmax = max([xmax varargin{i}.(xparam)]);
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Get the index of the nearest bin
for i=1:Ndata
  xidmin(i, 1) = nearest(varargin{i}.(xparam), xmin);
  xidmax(i, 1) = nearest(varargin{i}.(xparam), xmax);
end

if strcmp('freq',yparam) && strcmp('freq',dtype)
  tmpcfg = keepfields(cfg, {'parameter'});
  tmpcfg.avgoverfreq = 'yes';
  tmpcfg.frequency   = cfg.frequency;
  [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
  % restore the provenance information
  [cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
elseif strcmp('time',yparam) && strcmp('freq',dtype)
  tmpcfg = keepfields(cfg, {'parameter'});
  tmpcfg.avgovertime = 'yes';
  tmpcfg.latency     = cfg.latency;
  [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
  % restore the provenance information
  [cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
end

% Get physical y-axis range (ylim / parameter):
if strcmp(cfg.ylim, 'maxmin') || strcmp(cfg.ylim, 'maxabs')
  % Find maxmin throughout all varargins:
  ymin = [];
  ymax = [];
  for i=1:length(varargin)
    % Select the channels in the data that match with the layout and that
    % are selected for plotting:
    dat = [];
    dat = varargin{i}.(cfg.parameter);
    seldat1 = match_str(varargin{i}.label, lay.label);   % indexes labels corresponding in input and layout
    seldat2 = match_str(varargin{i}.label, cfg.channel); % indexes labels corresponding in input and plot-selection
    if isempty(seldat1)
      error('labels in data and labels in layout do not match');
    end
    data = dat(intersect(seldat1, seldat2), :);
    ymin = min([ymin min(min(min(data)))]);
    ymax = max([ymax max(max(max(data)))]);
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

% convert the layout to Ole's style of variable names
X      = lay.pos(:, 1);
Y      = lay.pos(:, 2);
width  = lay.width;
height = lay.height;
Lbl    = lay.label;

% Create empty channel coordinates and labels arrays:
chanX(1:length(Lbl)) = NaN;
chanY(1:length(Lbl)) = NaN;
chanLabels = cell(1, length(Lbl));

hold on;
colorLabels = [];

% Plot each data set:
for i=1:Ndata
  % Make vector dat with one value for each channel
  dat  = varargin{i}.(cfg.parameter);
  % get dimord dimensions
  dims = textscan(varargin{i}.dimord, '%s', 'Delimiter', '_');
  dims = dims{1};
  ydim = find(strcmp(yparam, dims));
  xdim = find(strcmp(xparam, dims));
  zdim = setdiff(1:ndims(dat), [ydim xdim]);
  % and permute
  dat = permute(dat, [zdim(:)' ydim xdim]);

  xval = varargin{i}.(xparam);

  % Take subselection of channels, this only works
  % in the non-interactive mode
  if exist('selchannel', 'var')
    sellab = match_str(varargin{i}.label, selchannel);
    label  = varargin{i}.label(sellab);
  else
    sellab = 1:numel(varargin{i}.label);
    label  = varargin{i}.label;
  end

  if isfull
    dat = dat(sel1, sel2, xidmin(i):xidmax(i));
    dat = nanmean(dat, meandir);
  elseif haslabelcmb
    dat = dat(sellab, xidmin(i):xidmax(i));
  else
    dat = dat(sellab, xidmin(i):xidmax(i));
  end
  xval = xval(xidmin(i):xidmax(i));

  % Select the channels in the data that match with the layout:
  [seldat, sellay] = match_str(label, cfg.layout.label);
  if isempty(seldat)
    error('labels in data and labels in layout do not match');
  end

  % gather the data of multiple input arguments
  datamatrix{i} = dat(seldat, :);

  % Select x and y coordinates and labels of the channels in the data
  layX = cfg.layout.pos(sellay, 1);
  layY = cfg.layout.pos(sellay, 2);
  layLabels = cfg.layout.label(sellay);

  if ~isempty(cfg.maskparameter)
    % one value for each channel, or one value for each channel-time point
    maskmatrix = varargin{1}.(cfg.maskparameter)(seldat, :);
    maskmatrix = maskmatrix(:, xidmin:xidmax);
  else
    % create an Nx0 matrix
    maskmatrix = zeros(length(seldat), 0);
  end

  if Ndata > 1
    if ischar(GRAPHCOLOR);        colorLabels = [colorLabels iname{i+1} '=' GRAPHCOLOR(i+1) '\n'];
    elseif isnumeric(GRAPHCOLOR); colorLabels = [colorLabels iname{i+1} '=' num2str(GRAPHCOLOR(i+1, :)) '\n'];
    end
  end
end % for number of input data

for m=1:length(layLabels)
  % Plot ER

  if ischar(GRAPHCOLOR);        color = GRAPHCOLOR(2:end);
  elseif isnumeric(GRAPHCOLOR); color = GRAPHCOLOR(2:end, :);
  end

  mask = maskmatrix(m, :);

  for i=1:Ndata
    yval(i, :) = datamatrix{i}(m, :);
  end

  % Clip out of bounds y values:
  yval(yval > ymax) = ymax;
  yval(yval < ymin) = ymin;

  if strcmp(cfg.showlabels, 'yes')
    label = layLabels(m);
  else
    % don't show labels
    label = [];
  end

  ft_plot_vector(xval, yval, 'width', width(m), 'height', height(m), 'hpos', layX(m), 'vpos', layY(m), 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'color', color, 'style', cfg.linestyle{i}, 'linewidth', cfg.linewidth, 'axis', cfg.axes, 'highlight', mask, 'highlightstyle', cfg.maskstyle, 'label', label, 'box', cfg.box, 'fontsize', cfg.fontsize);

  if i==1,
    % Keep ER plot coordinates (at centre of ER plot), and channel labels (will be stored in the figure's UserData struct):
    chanX(m) = X(m) + 0.5 * width(m);
    chanY(m) = Y(m) + 0.5 * height(m);
    chanLabels{m} = Lbl{m};
  end
end % for number of channels


% Add the colors of the different datasets to the comment:
cfg.comment = [cfg.comment colorLabels];

% Write comment text:
l = cellstrmatch('COMNT', Lbl);
if ~isempty(l)
  ft_plot_text(X(l), Y(l), sprintf(cfg.comment), 'Fontsize', cfg.fontsize, 'interpreter', 'none');
end

% Plot scales:
l = cellstrmatch('SCALE', Lbl);
if ~isempty(l)
  plotScales([xmin xmax], [ymin ymax], X(l), Y(l), width(1), height(1), cfg)
end

% set the figure window title
if isempty(get(gcf, 'Name'))
  if nargin > 1
    dataname = {inputname(2)};
    for k = 2:Ndata
      dataname{end+1} = inputname(k+1);
    end
  else % data provided through cfg.inputfile
    dataname = cfg.inputfile;
  end

  if isempty(cfg.figurename)
    set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), mfilename, join_str(', ', dataname)));
    set(gcf, 'NumberTitle', 'off');
  else
    set(gcf, 'name', cfg.figurename);
    set(gcf, 'NumberTitle', 'off');
  end
else
  dataname = {};
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')

  % add the dataname and channel information to the figure
  % this is used in the callbacks
  info          = guidata(gcf);
  info.x        = lay.pos(:, 1);
  info.y        = lay.pos(:, 2);
  info.label    = lay.label;
  info.dataname = dataname;
  guidata(gcf, info);

  set(gcf, 'WindowButtonUpFcn',  {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonMotionFcn'});
end

axis tight
axis off
if strcmp(cfg.box, 'yes')
  abc = axis;
  axis(abc + [-1 +1 -1 +1]*mean(abs(abc))/10)
end
hold off

% Set orientation for printing if specified
if ~isempty(cfg.orient)
  orient(gcf, cfg.orient);
end

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance

% add a menu to the figure, but only if the current figure does not have subplots
% also, delete any possibly existing previous menu, this is safe because delete([]) does nothing
delete(findobj(gcf, 'type', 'uimenu', 'label', 'FieldTrip'));
if numel(findobj(gcf, 'type', 'axes')) <= 1
  ftmenu = uimenu(gcf, 'Label', 'FieldTrip');
  uimenu(ftmenu, 'Label', 'Show pipeline', 'Callback', {@menu_pipeline, cfg});
  uimenu(ftmenu, 'Label', 'About', 'Callback', @menu_about);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotScales(hlim, vlim, hpos, vpos, width, height, cfg)

% the placement of all elements is identical
placement = {'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim};

ft_plot_box([hlim vlim], placement{:}, 'edgecolor', 'k');

if hlim(1)<=0 && hlim(2)>=0
  ft_plot_vector([0 0], vlim, placement{:}, 'color', 'b');
end

if vlim(1)<=0 && vlim(2)>=0
  ft_plot_vector(hlim, [0 0], placement{:}, 'color', 'b');
end

ft_plot_text(hlim(1), vlim(1), [num2str(hlim(1), 3) ' '], placement{:}, 'rotation', 90, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'top', 'Fontsize', cfg.fontsize);
ft_plot_text(hlim(2), vlim(1), [num2str(hlim(2), 3) ' '], placement{:}, 'rotation', 90, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'top', 'Fontsize', cfg.fontsize);
ft_plot_text(hlim(1), vlim(1), [num2str(vlim(1), 3) ' '], placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom', 'Fontsize', cfg.fontsize);
ft_plot_text(hlim(1), vlim(2), [num2str(vlim(2), 3) ' '], placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom', 'Fontsize', cfg.fontsize);

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
% SUBFUNCTION which is called after selecting channels in case of cfg.refchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_multiplotER(label, cfg, varargin)
if ~isempty(label)
  if isfield(cfg, 'inputfile')
    % the reading has already been done and varargin contains the data
    cfg = rmfield(cfg, 'inputfile');
  end
  % put data name in here, this cannot be resolved by other means
  info = guidata(gcf);
  cfg.dataname = info.dataname;
  if iscell(label)
    label = label{1};
  end
  cfg.refchannel = label; % FIXME this only works with label being a string
  fprintf('selected cfg.refchannel = ''%s''\n', cfg.refchannel);
  p = get(gcf, 'Position');
  f = figure;
  set(f, 'Position', p);
  ft_multiplotER(cfg, varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label, cfg, varargin)
if ~isempty(label)
  if isfield(cfg, 'inputfile')
    % the reading has already been done and varargin contains the data
    cfg = rmfield(cfg, 'inputfile');
  end
  cfg.channel = label;
  % put data name in here, this cannot be resolved by other means
  info = guidata(gcf);
  cfg.dataname = info.dataname;
  fprintf('selected cfg.channel = {');
  for i=1:(length(cfg.channel)-1)
    fprintf('''%s'', ', cfg.channel{i});
  end
  fprintf('''%s''}\n', cfg.channel{end});
  p = get(gcf, 'Position');
  f = figure;
  set(f, 'Position', p);
  ft_singleplotER(cfg, varargin{:});
end
