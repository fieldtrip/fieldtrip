function cfg = topoplot_common(cfg, varargin)

% TOPOPLOT_COMMON is shared by FT_TOPOPLOTTFR, FT_TOPOPLOTER and FT_TOPOPLOTIC, which
% serve as placeholder for the documentation and for the pre/postamble.

% Copyright (C) 2005-2011, F.C. Donders Centre
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'unused',     {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'cohrefchannel' 'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'deprecated', {'xparam'});

Ndata = numel(varargin);
if isnumeric(varargin{end})
  % the call with multiple inputs is done by ft_topoplotIC and recursively by ft_topoplotER/TFR itself
  % the last input argument points to the specific data to be plotted
  Ndata = Ndata - 1;
  indx  = varargin{end};
elseif isstruct(varargin{end})
  indx  = 1;
end

% the call with multiple inputs is done by ft_topoplotIC and recursively by ft_topoplotER/TFR itself
% the last input argument points to the specific data to be plotted
if Ndata>1 && ~isnumeric(varargin{end})
  for k=1:Ndata
    
    if k>1
      % create a new figure for the additional input arguments
      % ensure that the new figure appears at the same position
      figure('Position', get(gcf, 'Position'), 'Visible', get(gcf, 'Visible'));
    end
    
    % the indexing is necessary if ft_topoplotER/TFR is called from
    % ft_singleplotER/TFR when more input data structures exist. Somehow we need to
    % keep track of which of the data arguments is to be plotted (otherwise the first
    % data argument is only plotted). Yet, we cannot throw away the other data
    % structures, because in the interactive mode ft_singleplotER/TFR needs all data
    % again and the entry into ft_singleplotER/TFR will be through one of the figures,
    % which thus needs to have all data avalaible. At the moment I couldn't think of
    % anything better than using an additional indx variable and letting the function
    % recursively call itself.
    if isfield(cfg, 'inputfile')
      tmpcfg = rmfield(cfg, 'inputfile');
    else
      tmpcfg = cfg;
    end
    topoplot_common(tmpcfg, varargin{1:Ndata}, indx);
    indx = indx + 1;
  end
  return
end

if iscell(cfg.dataname)
  dataname = cfg.dataname{indx};
else
  dataname = cfg.dataname;
end

data = varargin{indx};
data = ft_checkdata(data, 'datatype', {'comp', 'timelock', 'freq'});

% check for option-values to be renamed
cfg = ft_checkconfig(cfg, 'renamedval', {'electrodes',     'dotnum',      'numbers'});
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim',           'absmax',      'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedback',    'inflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'highlight',      'yes',         'on'});

% check for renamed options
cfg = ft_checkconfig(cfg, 'renamed',     {'matrixside',    'directionality'});
cfg = ft_checkconfig(cfg, 'renamed',     {'electrodes',    'marker'});
cfg = ft_checkconfig(cfg, 'renamed',     {'emarker',       'markersymbol'});
cfg = ft_checkconfig(cfg, 'renamed',     {'ecolor',        'markercolor'});
cfg = ft_checkconfig(cfg, 'renamed',     {'emarkersize',   'markersize'});
cfg = ft_checkconfig(cfg, 'renamed',     {'efontsize',     'markerfontsize'});
cfg = ft_checkconfig(cfg, 'renamed',     {'hlmarker',      'highlightsymbol'});
cfg = ft_checkconfig(cfg, 'renamed',     {'hlcolor',       'highlightcolor'});
cfg = ft_checkconfig(cfg, 'renamed',     {'hlmarkersize',  'highlightsize'});
cfg = ft_checkconfig(cfg, 'renamed',     {'maplimits',     'zlim'});
% old ft_checkconfig adapted partially from topoplot.m (backwards backwards compatability)
cfg = ft_checkconfig(cfg, 'renamed',     {'grid_scale',    'gridscale'});
cfg = ft_checkconfig(cfg, 'renamed',     {'interpolate',   'interpolation'});
cfg = ft_checkconfig(cfg, 'renamed',     {'numcontour',    'contournum'});
cfg = ft_checkconfig(cfg, 'renamed',     {'electrod',      'marker'});
cfg = ft_checkconfig(cfg, 'renamed',     {'electcolor',    'markercolor'});
cfg = ft_checkconfig(cfg, 'renamed',     {'emsize',        'markersize'});
cfg = ft_checkconfig(cfg, 'renamed',     {'efsize',        'markerfontsize'});
cfg = ft_checkconfig(cfg, 'renamed',     {'headlimits',    'interplimits'});
% check for forbidden options
cfg = ft_checkconfig(cfg, 'forbidden',  {'hllinewidth', ...
  'headcolor', ...
  'hcolor', ...
  'hlinewidth', ...
  'contcolor', ...
  'outline', ...
  'highlightfacecolor', ...
  'showlabels'});

if ft_platform_supports('griddata-v4')
  default_interpmethod = 'v4';
else
  % Octave does not support 'v4', and 'cubic' not yet implemented
  default_interpmethod = 'linear';
end

% Set other config defaults
cfg.xlim              = ft_getopt(cfg, 'xlim',             'maxmin');
cfg.ylim              = ft_getopt(cfg, 'ylim',             'maxmin');
cfg.zlim              = ft_getopt(cfg, 'zlim',             'maxmin');
cfg.style             = ft_getopt(cfg, 'style',            'both');
cfg.gridscale         = ft_getopt(cfg, 'gridscale',         67);
cfg.interplimits      = ft_getopt(cfg, 'interplimits',     'head');
cfg.interpolation     = ft_getopt(cfg, 'interpolation',     default_interpmethod);
cfg.contournum        = ft_getopt(cfg, 'contournum',        6);
cfg.colorbar          = ft_getopt(cfg, 'colorbar',         'no');
cfg.shading           = ft_getopt(cfg, 'shading',          'flat');
cfg.comment           = ft_getopt(cfg, 'comment',          'auto');
cfg.commentpos        = ft_getopt(cfg, 'commentpos',        []);  % default is handled further down
cfg.fontsize          = ft_getopt(cfg, 'fontsize',          8);
cfg.fontweight        = ft_getopt(cfg, 'fontweight',       'normal');
cfg.baseline          = ft_getopt(cfg, 'baseline',         'no'); % to avoid warning in timelock/freqbaseline
cfg.trials            = ft_getopt(cfg, 'trials',           'all', 1);
cfg.interactive       = ft_getopt(cfg, 'interactive',      'yes');
cfg.hotkeys           = ft_getopt(cfg, 'hotkeys',          'yes');
cfg.renderer          = ft_getopt(cfg, 'renderer',          []); % MATLAB sets the default
cfg.marker            = ft_getopt(cfg, 'marker',           'on');
cfg.markersymbol      = ft_getopt(cfg, 'markersymbol',     'o');
cfg.markercolor       = ft_getopt(cfg, 'markercolor',       [0 0 0]);
cfg.markersize        = ft_getopt(cfg, 'markersize',        2);
cfg.markerfontsize    = ft_getopt(cfg, 'markerfontsize',    8);
cfg.highlight         = ft_getopt(cfg, 'highlight',        'off');
cfg.highlightchannel  = ft_getopt(cfg, 'highlightchannel', 'all', 1); % highlight may be 'on', making highlightchannel {} meaningful
cfg.highlightsymbol   = ft_getopt(cfg, 'highlightsymbol',  '*');
cfg.highlightcolor    = ft_getopt(cfg, 'highlightcolor',    [0 0 0]);
cfg.highlightsize     = ft_getopt(cfg, 'highlightsize',     6);
cfg.highlightfontsize = ft_getopt(cfg, 'highlightfontsize', 8);
cfg.labeloffset       = ft_getopt(cfg, 'labeloffset',       0.005);
cfg.maskparameter     = ft_getopt(cfg, 'maskparameter',     []);
cfg.component         = ft_getopt(cfg, 'component',         []);
cfg.directionality    = ft_getopt(cfg, 'directionality',    []);
cfg.channel           = ft_getopt(cfg, 'channel',          'all');
cfg.refchannel        = ft_getopt(cfg, 'refchannel',        []);
cfg.figurename        = ft_getopt(cfg, 'figurename',        []);
cfg.interpolatenan    = ft_getopt(cfg, 'interpolatenan',   'yes');

% default commentpos
if isempty(cfg.commentpos)
  if any(ismember(cfg.layout.label, 'COMNT'))
    % layout went through ft_prepare_layout in ft_topoplotER/TFR
    % use the position that is specified in the layout
    cfg.commentpos = 'layout';
  else
    % put it in the left bottom
    cfg.commentpos = 'leftbottom';
  end
end

% the user can either specify a single group of channels for highlighting (all in the
% same style), or multiple groups with a different style for each group. The latter
% is used by ft_clusterplot.
if ~iscell(cfg.highlight)
  cfg.highlight = {cfg.highlight};
end
if iscell(cfg.highlightchannel) && ~(iscell(cfg.highlightchannel{1}) || isnumeric(cfg.highlightchannel{1}))
  cfg.highlightchannel = {cfg.highlightchannel};
elseif ischar(cfg.highlightchannel)
  cfg.highlightchannel = {{cfg.highlightchannel}};
end
if ~iscell(cfg.highlightsymbol)
  cfg.highlightsymbol = {cfg.highlightsymbol};
end
if ~iscell(cfg.highlightcolor)
  cfg.highlightcolor = {cfg.highlightcolor};
end
if ~iscell(cfg.highlightsize)
  cfg.highlightsize = {cfg.highlightsize};
end
if ~iscell(cfg.highlightfontsize)
  cfg.highlightfontsize = {cfg.highlightfontsize};
end
% then make sure all cell-arrays for options have length ncellhigh and default the last element if not present
ncellhigh = length(cfg.highlightchannel);
if length(cfg.highlightsymbol)    < ncellhigh,   cfg.highlightsymbol{ncellhigh}    = 'o';       end
if length(cfg.highlightcolor)     < ncellhigh,   cfg.highlightcolor{ncellhigh}     = [0 0 0];   end
if length(cfg.highlightsize)      < ncellhigh,   cfg.highlightsize{ncellhigh}      = 6;         end
if length(cfg.highlightfontsize)  < ncellhigh,   cfg.highlightfontsize{ncellhigh}  = 8;         end
% then default all empty cells
for icell = 1:ncellhigh
  if isempty(cfg.highlightsymbol{icell}),    cfg.highlightsymbol{icell} = 'o';     end
  if isempty(cfg.highlightcolor{icell}),     cfg.highlightcolor{icell} = [0 0 0];  end
  if isempty(cfg.highlightsize{icell}),      cfg.highlightsize{icell} = 6;         end
  if isempty(cfg.highlightfontsize{icell}),  cfg.highlightfontsize{icell} = 8;     end
end

% for backwards compatability
if strcmp(cfg.marker, 'highlights')
  ft_warning('using cfg.marker option -highlights- is no longer used, please use cfg.highlight')
  cfg.marker = 'off';
end

% check colormap is proper format and set it
if isfield(cfg, 'colormap')
  if ~isnumeric(cfg.colormap)
    cfg.colormap = colormap(cfg.colormap);
  end
  if size(cfg.colormap,2)~=3
    ft_error('cfg.colormap must be Nx3');
  end
  colormap(cfg.colormap);
  ncolors = size(cfg.colormap,1);
else
  ncolors = []; % let the low-level function deal with this
end


%% Section 2: data handling, this also includes converting bivariate (chan_chan and chancmb) into univariate data

dtype  = ft_datatype(data);
hastime = isfield(data, 'time');

% Set x/y/parameter defaults according to datatype and dimord
switch dtype
  case 'timelock'
    xparam = 'time';
    yparam = '';
    if isfield(data, 'trial')
      cfg.parameter = ft_getopt(cfg, 'parameter', 'trial');
    elseif isfield(data, 'individual')
      cfg.parameter = ft_getopt(cfg, 'parameter', 'individual');
    elseif isfield(data, 'avg')
      cfg.parameter = ft_getopt(cfg, 'parameter', 'avg');
    end
  case 'freq'
    if hastime
      xparam = 'time';
      yparam = 'freq';
      cfg.parameter = ft_getopt(cfg, 'parameter', 'powspctrm');
    else
      xparam = 'freq';
      yparam = '';
      cfg.parameter = ft_getopt(cfg, 'parameter', 'powspctrm');
    end
  case 'comp'
    % Add a pseudo-axis with the component numbers
    data.comp = 1:size(data.topo,2);
    if ~isempty(cfg.component)
      % make a selection of components
      data.comp  = data.comp(cfg.component);
      data.topo  = data.topo(:,cfg.component);
      try, data.label     = data.label(cfg.component); end
      try, data.unmixing  = data.unmixing(cfg.component,:); end
    end
    % Rename the field with topographic label information
    data.label      = data.topolabel;
    data.topodimord = 'chan_comp';
    data = removefields(data, {'topolabel', 'unmixing', 'unmixingdimord'}); % not needed any more
    xparam = 'comp';
    yparam = '';
    cfg.parameter = ft_getopt(cfg, 'parameter', 'topo');
  otherwise
    % if the input data is not one of the standard data types, or if the functional
    % data is just one value per channel: in this case xparam, yparam are not defined
    % and the user should define the parameter
    if ~isfield(data, 'label'),     ft_error('the input data should at least contain a label-field');            end
    if ~isfield(cfg,  'parameter'), ft_error('the configuration should at least contain a ''parameter'' field'); end
    if ~isfield(cfg,  'xparam')
      cfg.xlim = [1 1];
      xparam   = '';
    end
end

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

% Apply baseline correction
if ~strcmp(cfg.baseline, 'no')
  % keep mask-parameter if it is set
  if ~isempty(cfg.maskparameter)
    tempmask = data.(cfg.maskparameter);
  end
  if strcmp(xparam, 'time') && strcmp(yparam, 'freq')
    data = ft_freqbaseline(cfg, data);
  elseif strcmp(xparam, 'time') && strcmp(yparam, '')
    data = ft_timelockbaseline(cfg, data);
  end
  % put mask-parameter back if it is set
  if ~isempty(cfg.maskparameter)
    data.(cfg.maskparameter) = tempmask;
  end
end

% time and/or frequency should NOT be selected and averaged here, since a singleplot might follow in interactive mode
tmpcfg = keepfields(cfg, {'channel', 'showcallinfo', 'trials'});
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
dimord = getdimord(varargin{1}, cfg.parameter);
if startsWith(dimord, 'chan_chan_') || startsWith(dimord, 'chancmb_')
  % convert the bivariate data to univariate and call the parent plotting function again
  s = dbstack;
  cfg.originalfunction = s(2).name;
  cfg.trials = 'all'; % trial selection has been taken care off
  bivariate_common(cfg, varargin{:});
  return
end

% Apply channel-type specific scaling
tmpcfg = keepfields(cfg, {'parameter', 'chanscale', 'ecgscale', 'eegscale', 'emgscale', 'eogscale', 'gradscale', 'magscale', 'megscale', 'mychan', 'mychanscale'});
data = chanscale_common(tmpcfg, data);


%% Section 3: select the data to be plotted and determine min/max range

dimord = getdimord(varargin{1}, cfg.parameter);
dimtok = tokenize(dimord, '_');

% Create time-series of small topoplots
if ~ischar(cfg.xlim) && length(cfg.xlim)>2 %&& any(ismember(dimtok, 'time'))
  % Switch off interactive mode:
  cfg.interactive = 'no';
  xlims = cfg.xlim;
  % Iteratively call topoplotER with different xlim values:
  nplots = numel(xlims)-1;
  nyplot = ceil(sqrt(nplots));
  nxplot = ceil(nplots./nyplot);
  for i=1:length(xlims)-1
    subplot(nxplot, nyplot, i);
    cfg.xlim = xlims(i:i+1);
    ft_topoplotTFR(cfg, data);
  end
  return
end

% Get physical min/max range of x
if strcmp(cfg.xlim, 'maxmin')
  xmin = min(data.(xparam));
  xmax = max(data.(xparam));
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end
xminindx = nearest(data.(xparam), xmin);
xmaxindx = nearest(data.(xparam), xmax);
xmin = data.(xparam)(xminindx);
xmax = data.(xparam)(xmaxindx);
selx = xminindx:xmaxindx;

% Get physical min/max range of y
if ~isempty(yparam)
  if strcmp(cfg.ylim, 'maxmin')
    ymin = min(data.(yparam));
    ymax = max(data.(yparam));
  else
    ymin = cfg.ylim(1);
    ymax = cfg.ylim(2);
  end
  yminindx = nearest(data.(yparam), ymin);
  ymaxindx = nearest(data.(yparam), ymax);
  ymin = data.(yparam)(yminindx);
  ymax = data.(yparam)(ymaxindx);
  sely = yminindx:ymaxindx;
end

% Take subselection of channels, this only works if the interactive mode is switched off
if exist('selchannel', 'var')
  sellab = match_str(data.label, selchannel);
  label  = data.label(sellab);
else
  sellab = 1:numel(data.label);
  label  = data.label;
end

% Make data vector with one scalar value for each channel
dat = data.(cfg.parameter);
% get dimord dimensions
ydim = find(strcmp(yparam, dimtok));
xdim = find(strcmp(xparam, dimtok));
zdim = setdiff(1:ndims(dat), [ydim xdim]);
% and permute
dat = permute(dat, [zdim(:)' ydim xdim]);

if ~isempty(yparam)
  % time-frequency data
  dat = dat(sellab, sely, selx);
  dat = nanmean(nanmean(dat, 3), 2);
elseif ~isempty(cfg.component)
  % component data, nothing to do
else
  % time or frequency data
  dat = dat(sellab, selx);
  dat = nanmean(dat, 2);
end
dat = dat(:);

if isfield(data, cfg.maskparameter)
  % Make mask vector with one value for each channel
  msk = data.(cfg.maskparameter);
  % get dimord dimensions
  ydim = find(strcmp(yparam, dimtok));
  xdim = find(strcmp(xparam, dimtok));
  zdim = setdiff(1:ndims(dat), [ydim xdim]);
  % and permute
  msk = permute(msk, [zdim(:)' ydim xdim]);
  
  if ~isempty(yparam)
    % time-frequency data
    msk = msk(sellab, sely, selx);
  elseif ~isempty(cfg.component)
    % component data, nothing to do
  else
    % time or frequency data
    msk = msk(sellab, selx);
  end
  
  if size(msk,2)>1 || size(msk,3)>1
    ft_warning('no masking possible for average over multiple latencies or frequencies -> cfg.maskparameter cleared')
    msk = [];
  end
  
else
  msk = [];
end

% Select the channels in the data that match with the layout:
[seldat, sellay] = match_str(label, cfg.layout.label);
if isempty(seldat)
  ft_error('labels in data and labels in layout do not match');
end

dat = dat(seldat);
if ~isempty(msk)
  msk = msk(seldat);
end

% Select x and y coordinates and labels of the channels in the data
chanX = cfg.layout.pos(sellay,1);
chanY = cfg.layout.pos(sellay,2);

% Get physical min/max range of z:
if strcmp(cfg.zlim, 'maxmin')
  zmin = min(dat);
  zmax = max(dat);
elseif strcmp(cfg.zlim, 'maxabs')
  zmin = -max(max(abs(dat)));
  zmax = max(max(abs(dat)));
elseif strcmp(cfg.zlim, 'zeromax')
  zmin = 0;
  zmax = max(dat);
elseif strcmp(cfg.zlim, 'minzero')
  zmin = min(dat);
  zmax = 0;
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% make comment
switch cfg.comment
  case 'no'
    cfg.comment = '';
  case 'auto'
    comment = date;
    if ~isempty(xparam)
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, xparam, xmin, xmax);
    end
    if ~isempty(yparam)
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, yparam, ymin, ymax);
    end
    if ~isempty(cfg.parameter)
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, cfg.parameter, zmin, zmax);
    end
    cfg.comment = comment;
  case 'xlim'
    comment = date;
    if ~isempty(xparam)
      cfg.comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, xparam, xmin, xmax);
    end
  case 'ylim'
    comment = date;
    if ~isempty(yparam)
      cfg.comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, yparam, ymin, ymax);
    end
  case 'zlim'
    comment = date;
    if ~isempty(yparam)
      cfg.comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, cfg.parameter, zmin, zmax);
    end
end % switch comment

if ~isempty(cfg.refchannel)
  if iscell(cfg.refchannel)
    cfg.comment = sprintf('%s\nreference=%s %s', cfg.comment, cfg.refchannel{:});
  else
    cfg.comment = sprintf('%s\nreference=%s %s', cfg.comment, cfg.refchannel);
  end
end

% Specify the x and y coordinates of the comment
if strcmp(cfg.commentpos, 'layout')
  ind_comment = find(strcmp(cfg.layout.label, 'COMNT'));
  x_comment = cfg.layout.pos(ind_comment,1);
  y_comment = cfg.layout.pos(ind_comment,2);
elseif strcmp(cfg.commentpos, 'lefttop')
  x_comment = -0.7;
  y_comment =  0.6;
  HorAlign = 'left';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos, 'leftbottom')
  x_comment = -0.6;
  y_comment = -0.6;
  HorAlign = 'left';
  VerAlign = 'bottom';
elseif strcmp(cfg.commentpos, 'middletop')
  x_comment =  0;
  y_comment =  0.75;
  HorAlign = 'center';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos, 'middlebottom')
  x_comment =  0;
  y_comment = -0.7;
  HorAlign = 'center';
  VerAlign = 'bottom';
elseif strcmp(cfg.commentpos, 'righttop')
  x_comment =  0.65;
  y_comment =  0.6;
  HorAlign = 'right';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos, 'rightbottom')
  x_comment =  0.6;
  y_comment = -0.6;
  HorAlign = 'right';
  VerAlign = 'bottom';
elseif isnumeric(cfg.commentpos)
  x_comment = cfg.commentpos(1);
  y_comment = cfg.commentpos(2);
  HorAlign = 'left';
  VerAlign = 'middle';
  x_comment = 0.9*((x_comment-min(cfg.xlim))/(max(cfg.xlim)-min(cfg.xlim))-0.5);
  y_comment = 0.9*((y_comment-min(cfg.ylim))/(max(cfg.ylim)-min(cfg.ylim))-0.5);
end

% Draw topoplot
cla
hold on

% check for nans
nanInds = isnan(dat);
if strcmp(cfg.interpolatenan, 'yes') && any(nanInds)
  ft_warning('removing NaNs from the data');
  chanX(nanInds) = [];
  chanY(nanInds) = [];
  dat(nanInds)   = [];
  if ~isempty(msk)
    msk(nanInds) = [];
  end
end

% Set ft_plot_topo specific options
if strcmp(cfg.interplimits, 'head')
  interplimits = 'mask';
else
  interplimits = cfg.interplimits;
end
if strcmp(cfg.style, 'both');            style = 'surfiso';     end
if strcmp(cfg.style, 'straight');        style = 'surf';        end
if strcmp(cfg.style, 'contour');         style = 'iso';         end
if strcmp(cfg.style, 'fill');            style = 'isofill';     end
if strcmp(cfg.style, 'straight_imsat');  style = 'imsat';       end
if strcmp(cfg.style, 'both_imsat');      style = 'imsatiso';    end

% Draw plot
if ~strcmp(cfg.style, 'blank')
  opt = {'interpmethod', cfg.interpolation, ...
    'interplim',    interplimits, ...
    'gridscale',    cfg.gridscale, ...
    'outline',      cfg.layout.outline, ...
    'shading',      cfg.shading, ...
    'isolines',     cfg.contournum, ...
    'mask',         cfg.layout.mask, ...
    'style',        style, ...
    'datmask',      msk};
  if strcmp(style, 'imsat') || strcmp(style, 'imsatiso')
    % add clim to opt
    opt = [opt {'clim', [zmin zmax], 'ncolors',ncolors}];
  end
  ft_plot_topo(chanX, chanY, dat, opt{:});
elseif ~strcmp(cfg.style, 'blank')
  ft_plot_lay(cfg.layout, 'box', 'no', 'label', 'no', 'point', 'no')
end

% For Highlight (channel-selection)
for icell = 1:length(cfg.highlight)
  if ~strcmp(cfg.highlight{icell}, 'off')
    cfg.highlightchannel{icell} = ft_channelselection(cfg.highlightchannel{icell}, data.label);
    [dum, layoutindex] = match_str(cfg.highlightchannel{icell}, cfg.layout.label);
    templay = [];
    templay.outline = cfg.layout.outline;
    templay.mask    = cfg.layout.mask;
    templay.pos     = cfg.layout.pos(layoutindex,:);
    templay.width   = cfg.layout.width(layoutindex);
    templay.height  = cfg.layout.height(layoutindex);
    templay.label   = cfg.layout.label(layoutindex);
    if strcmp(cfg.highlight{icell}, 'labels') || strcmp(cfg.highlight{icell}, 'numbers')
      labelflg = 1;
    else
      labelflg = 0;
    end
    if strcmp(cfg.highlight{icell}, 'numbers')
      for ichan = 1:length(layoutindex)
        templay.label{ichan} = num2str(match_str(data.label, templay.label{ichan}));
      end
    end
    
    ft_plot_lay(templay, 'box', 'no', 'label', labelflg, 'point', ~labelflg, ...
      'pointsymbol',  cfg.highlightsymbol{icell}, ...
      'pointcolor',   cfg.highlightcolor{icell}, ...
      'pointsize',    cfg.highlightsize{icell}, ...
      'fontsize',     cfg.highlightfontsize{icell}, ...
      'labeloffset',  cfg.labeloffset, ...
      'labelalignh', 'center', ...
      'labelalignv', 'middle');
  end
end % for icell

% For Markers (all channels)
cfg = ft_checkopt(cfg, 'marker', {}, {'on', 'off', 'labels', 'numbers'});
if ~strcmp(cfg.marker, 'off')
  channelsToMark = 1:length(data.label);
  channelsToHighlight = [];
  for icell = 1:length(cfg.highlight)
    if ~strcmp(cfg.highlight{icell}, 'off')
      channelsToHighlight = [channelsToHighlight; match_str(data.label, cfg.highlightchannel{icell})];
    end
  end
  if strcmp(cfg.interpolatenan, 'no')
    channelsNotMark = channelsToHighlight;
  else
    channelsNotMark = union(find(isnan(dat)), channelsToHighlight);
  end
  channelsToMark(channelsNotMark) = [];
  [dum, layoutindex] = match_str(ft_channelselection(channelsToMark, data.label), cfg.layout.label);
  templay = [];
  templay.outline = cfg.layout.outline;
  templay.mask    = cfg.layout.mask;
  templay.pos     = cfg.layout.pos(layoutindex,:);
  templay.width   = cfg.layout.width(layoutindex);
  templay.height  = cfg.layout.height(layoutindex);
  templay.label   = cfg.layout.label(layoutindex);
  if strcmp(cfg.marker, 'labels') || strcmp(cfg.marker, 'numbers')
    labelflg = 1;
  else
    labelflg = 0;
  end
  if strcmp(cfg.marker, 'numbers')
    for ichan = 1:length(layoutindex)
      templay.label{ichan} = num2str(match_str(data.label,templay.label{ichan}));
    end
  end
  ft_plot_lay(templay, 'box', 'no', 'label',labelflg, 'point', ~labelflg, ...
    'pointsymbol',  cfg.markersymbol, ...
    'pointcolor',   cfg.markercolor, ...
    'pointsize',    cfg.markersize, ...
    'fontsize',     cfg.markerfontsize, ...
    'labeloffset',  cfg.labeloffset, ...
    'labelalignh', 'center', ...
    'labelalignv', 'middle');
end

if isfield(cfg, 'vector')
  % FIXME this is not documented
  vecX = nanmean(real(data.(cfg.vector)(:,selx)), 2);
  vecY = nanmean(imag(data.(cfg.vector)(:,selx)), 2);
  
  % scale quiver relative to largest gradiometer sample
  k = 0.15/max([max(abs(real(data.(cfg.vector)(:)))) max(abs(imag(data.(cfg.vector)(:))))]);
  quiver(chanX, chanY, k*vecX, k*vecY, 0, 'red');
end

% Write comment
if ~strcmp(cfg.comment, 'no')
  if strcmp(cfg.commentpos, 'title')
    title(cfg.comment, 'FontSize', cfg.fontsize);
  else
    ft_plot_text(x_comment, y_comment, cfg.comment, 'FontSize', cfg.fontsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontWeight', cfg.fontweight);
  end
end

% Set colour axis
if ~strcmp(cfg.style, 'blank')
  caxis([zmin zmax]);
end

% Plot colorbar
if isfield(cfg, 'colorbar')
  if strcmp(cfg.colorbar, 'yes')
    colorbar;
  elseif ~strcmp(cfg.colorbar, 'no')
    colorbar('location', cfg.colorbar);
  end
end

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

% set the figure window title, but only if the user has not changed it
if isempty(get(gcf, 'Name'))
  if isfield(cfg, 'funcname')
    funcname = cfg.funcname;
  else
    funcname = mfilename;
  end
  if isempty(cfg.figurename)
    dataname_str = join_str(', ', dataname);
    set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), funcname, dataname_str));
    set(gcf, 'NumberTitle', 'off');
  else
    set(gcf, 'name', cfg.figurename);
    set(gcf, 'NumberTitle', 'off');
  end
end

axis off
hold off
axis equal

if strcmp('yes', cfg.hotkeys)
  %  Attach data and cfg to figure and attach a key listener to the figure
  set(gcf, 'KeyPressFcn', {@key_sub, zmin, zmax})
end

% Make the figure interactive
if strcmp(cfg.interactive, 'yes')
  % add the cfg/data/channel information to the figure under identifier linked to this axis
  ident                    = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
  set(gca, 'tag',ident);
  info                     = guidata(gcf);
  info.(ident).x           = cfg.layout.pos(:, 1);
  info.(ident).y           = cfg.layout.pos(:, 2);
  info.(ident).label       = cfg.layout.label;
  info.(ident).dataname    = dataname;
  info.(ident).cfg         = cfg;
  if ~isfield(info.(ident),'datvarargin')
    info.(ident).datvarargin = varargin(1:Ndata); % add all datasets to figure
  end
  info.(ident).datvarargin{indx} = data; % update current dataset (e.g. baselined, channel selection, etc)
  guidata(gcf, info);
  if any(strcmp(data.dimord, {'chan_time', 'chan_freq', 'subj_chan_time', 'rpt_chan_time', 'chan_chan_freq', 'chancmb_freq', 'rpt_chancmb_freq', 'subj_chancmb_freq'}))
    set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER}, 'event', 'WindowButtonMotionFcn'});
  elseif any(strcmp(data.dimord, {'chan_freq_time', 'subj_chan_freq_time', 'rpt_chan_freq_time', 'rpttap_chan_freq_time', 'chan_chan_freq_time', 'chancmb_freq_time', 'rpt_chancmb_freq_time', 'subj_chancmb_freq_time'}))
    set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR}, 'event', 'WindowButtonMotionFcn'});
  else
    ft_warning('unsupported dimord "%s" for interactive plotting', data.dimord);
  end
end

% add a menu to the figure, but only if the current figure does not have subplots
% also, delete any possibly existing previous menu, this is safe because delete([]) does nothing
delete(findobj(gcf, 'type', 'uimenu', 'label', 'FieldTrip'));
if numel(findobj(gcf, 'type', 'axes')) <= 1
  ftmenu = uimenu(gcf, 'Label', 'FieldTrip');
  if ft_platform_supports('uimenu')
    % not supported by Octave
    uimenu(ftmenu, 'Label', 'Show pipeline',  'Callback', {@menu_pipeline, cfg});
    uimenu(ftmenu, 'Label', 'About',  'Callback', @menu_about);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label)
ident       = get(gca, 'tag');
info        = guidata(gcf);
cfg         = info.(ident).cfg;
datvarargin = info.(ident).datvarargin;
if ~isempty(label)
  cfg = removefields(cfg, 'inputfile');   % the reading has already been done and varargin contains the data
  cfg.baseline = 'no';                    % make sure the next function does not apply a baseline correction again
  cfg.channel = label;
  cfg.dataname = info.(ident).dataname;   % put data name in here, this cannot be resolved by other means
  cfg.trials = 'all';                     % trial selection has already been taken care of
  cfg.xlim = 'maxmin';
  % if user specified a zlim, copy it over to the ylim of singleplot
  if isfield(cfg, 'zlim')
    cfg.ylim = cfg.zlim;
    cfg = rmfield(cfg, 'zlim');
  end
  fprintf('selected cfg.channel = {%s}\n', join_str(', ', cfg.channel));
  % ensure that the new figure appears at the same position
  f = figure('Position', get(gcf, 'Position'), 'Visible', get(gcf, 'Visible'));
  ft_singleplotER(cfg, datvarargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotTFR(label)
ident       = get(gca, 'tag');
info        = guidata(gcf);
cfg         = info.(ident).cfg;
datvarargin = info.(ident).datvarargin;
if ~isempty(label)
  cfg = removefields(cfg, 'inputfile');   % the reading has already been done and varargin contains the data
  cfg.baseline = 'no';                    % make sure the next function does not apply a baseline correction again
  cfg.channel = label;
  cfg.dataname = info.(ident).dataname;   % put data name in here, this cannot be resolved by other means
  cfg.trials = 'all';                     % trial selection has already been taken care of
  cfg.xlim = 'maxmin';
  cfg.ylim = 'maxmin';
  fprintf('selected cfg.channel = {%s}\n', join_str(', ', cfg.channel));
  % ensure that the new figure appears at the same position
  f = figure('Position', get(gcf, 'Position'), 'Visible', get(gcf, 'Visible'));
  ft_singleplotTFR(cfg, datvarargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata, varargin)
climits = caxis;
incr_c  = abs(climits(2) - climits(1)) /10;

if length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:}, 'control')
  % TRANSLATE by 10%
  switch eventdata.Key
    case 'pageup'
      caxis([climits(1)+incr_c climits(2)+incr_c]);
    case 'pagedown'
      caxis([climits(1)-incr_c climits(2)-incr_c]);
  end % switch
else
  % ZOOM by 10%
  switch eventdata.Key
    case 'pageup'
      caxis([climits(1)-incr_c climits(2)+incr_c]);
    case 'pagedown'
      caxis([climits(1)+incr_c climits(2)-incr_c]);
    case 'm'
      caxis([varargin{1} varargin{2}]);
  end % switch
end % if
