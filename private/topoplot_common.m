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

revision = '$Id$';

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'unused',     {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'cohrefchannel' 'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'deprecated', {'xparam'});

Ndata = numel(varargin);
if ~isempty(varargin) && isnumeric(varargin{end})
  Ndata = Ndata - 1;
  indx  = varargin{end};
else
  indx  = 1;
end

% the call with multiple inputs is done by ft_topoplotIC and recursively by ft_topoplotTFR itself
if Ndata>1 && ~isnumeric(varargin{end})
  for k=1:Ndata
    
    if k>1
      % create a new figure for the additional input arguments
      % ensure that the new figure appears at the same position
      f = figure('Position', get(gcf, 'Position'), 'Visible', get(gcf, 'Visible'));
    end

    if isfield(cfg, 'inputfile')
      cfg = rmfield(cfg, 'inputfile');
    end
    
    % the indexing is necessary if ft_topoplotTFR is called from
    % ft_singleplotER when more input data structures exist. somehow we need to
    % keep track of which of the data arguments is to be plotted (otherwise the
    % first data argument is only plotted). yet, we cannot throw away the
    % other data structures, because in the interactive mode
    % ft_singleplotER needs all data again and the entry into
    % ft_singleplotER will be through one of the figures (which thus needs
    % to have all data avalaible. at the moment I couldn't think of
    % anything better than using an additional indx variable and letting the
    % function recursively call itself.
    topoplot_common(cfg, varargin{1:Ndata}, indx);
    indx = indx + 1;
  end
  return
end

data = varargin{indx};
data = ft_checkdata(data, 'datatype', {'comp', 'timelock', 'freq'});

% check for option-values to be renamed
cfg = ft_checkconfig(cfg, 'renamedval', {'electrodes',     'dotnum',      'numbers'});
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim',           'absmax',      'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedback',    'inflow'});

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
  default_interpmethod='v4';
else
  % Octave does not support 'v4', and 'cubic' not yet implemented
  default_interpmethod='linear';
end

% Set other config defaults
cfg.xlim              = ft_getopt(cfg, 'xlim',          'maxmin');
cfg.ylim              = ft_getopt(cfg, 'ylim',          'maxmin');
cfg.zlim              = ft_getopt(cfg, 'zlim',          'maxmin');
cfg.style             = ft_getopt(cfg, 'style',         'both');
cfg.gridscale         = ft_getopt(cfg, 'gridscale',     67);
cfg.interplimits      = ft_getopt(cfg, 'interplimits',  'head');
cfg.interpolation     = ft_getopt(cfg, 'interpolation', default_interpmethod);
cfg.contournum        = ft_getopt(cfg, 'contournum',    6);
cfg.colorbar          = ft_getopt(cfg, 'colorbar',      'no');
cfg.shading           = ft_getopt(cfg, 'shading',       'flat');
cfg.comment           = ft_getopt(cfg, 'comment',       'auto');
cfg.commentpos        = ft_getopt(cfg, 'commentpos',    []);  % default is handled further down
cfg.fontsize          = ft_getopt(cfg, 'fontsize',      8);
cfg.fontweight        = ft_getopt(cfg, 'fontweight',    'normal');
cfg.baseline          = ft_getopt(cfg, 'baseline',      'no'); % to avoid warning in timelock/freqbaseline
cfg.trials            = ft_getopt(cfg, 'trials',        'all', 1);
cfg.interactive       = ft_getopt(cfg, 'interactive',   'yes');
cfg.hotkeys           = ft_getopt(cfg, 'hotkeys',       'no');
cfg.renderer          = ft_getopt(cfg, 'renderer',      []); % MATLAB sets the default
cfg.marker            = ft_getopt(cfg, 'marker',        'on');
cfg.markersymbol      = ft_getopt(cfg, 'markersymbol',  'o');
cfg.markercolor       = ft_getopt(cfg, 'markercolor',   [0 0 0]);
cfg.markersize        = ft_getopt(cfg, 'markersize',    2);
cfg.markerfontsize    = ft_getopt(cfg, 'markerfontsize', 8);
cfg.highlight         = ft_getopt(cfg, 'highlight',     'off');
cfg.highlightchannel  = ft_getopt(cfg, 'highlightchannel',  'all', 1); % highlight may be 'on', making highlightchannel {} meaningful
cfg.highlightsymbol   = ft_getopt(cfg, 'highlightsymbol',   '*');
cfg.highlightcolor    = ft_getopt(cfg, 'highlightcolor',    [0 0 0]);
cfg.highlightsize     = ft_getopt(cfg, 'highlightsize',     6);
cfg.highlightfontsize = ft_getopt(cfg, 'highlightfontsize', 8);
cfg.labeloffset       = ft_getopt(cfg, 'labeloffset',       0.005);
cfg.maskparameter     = ft_getopt(cfg, 'maskparameter',     []);
cfg.component         = ft_getopt(cfg, 'component',         []);
cfg.directionality    = ft_getopt(cfg, 'directionality',    []);
cfg.channel           = ft_getopt(cfg, 'channel',           'all');
cfg.figurename        = ft_getopt(cfg, 'figurename',        []);
cfg.interpolatenan    = ft_getopt(cfg, 'interpolatenan',    'yes');

% default commentpos
if isempty(cfg.commentpos)
  if ~isempty(strcmp(cfg.layout.label, 'COMNT')) 
    % layout went through ft_prepare_layout in ft_topoplotER/TFR
    % use the position that is specified in the layout
    cfg.commentpos = 'layout';
  else
    % put it in the left bottom
    cfg.commentpos = 'leftbottom';
  end
end

% compatibility for previous highlighting option
if isnumeric(cfg.highlight)
  cfg.highlightchannel = cfg.highlight;
  cfg.highlight = 'on';
  warning('cfg.highlight is now used for specifying highlighting-mode, use cfg.highlightchannel instead of cfg.highlight for specifying channels')
elseif iscell(cfg.highlight)
  if ~iscell(cfg.highlightchannel)
    cfg.highlightchannel = cell(1,length(cfg.highlight));
  end
  for icell = 1:length(cfg.highlight)
    if isnumeric(cfg.highlight{icell})
      cfg.highlightchannel{icell} = cfg.highlight{icell};
      cfg.highlight{icell} = 'on';
      warning('cfg.highlight is now used for specifying highlighting-mode, use cfg.highlightchannel instead of cfg.highlight for specifying channels')
    end
  end
end

% Converting all highlight options to cell-arrays if they're not cell-arrays,
% to make defaulting, checking for backwards compatibility and error
% checking easier
if ~iscell(cfg.highlight),            cfg.highlight         = {cfg.highlight};            end
if isempty(cfg.highlightchannel),     cfg.highlightchannel  = ''; end
if ~iscell(cfg.highlightchannel),     cfg.highlightchannel  = {cfg.highlightchannel};     end
if ischar(cfg.highlightchannel{1}),   cfg.highlightchannel  = {cfg.highlightchannel};     end % {'all'} is valid input to channelselection, {1:5} isn't
if ~iscell(cfg.highlightsymbol),      cfg.highlightsymbol   = {cfg.highlightsymbol};      end
if ~iscell(cfg.highlightcolor),       cfg.highlightcolor    = {cfg.highlightcolor};       end
if ~iscell(cfg.highlightsize),        cfg.highlightsize     = {cfg.highlightsize};        end
if ~iscell(cfg.highlightfontsize),    cfg.highlightfontsize = {cfg.highlightfontsize};    end
% then make sure all cell-arrays for options have length ncellhigh and default the last element if not present
ncellhigh = length(cfg.highlight);
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
if strcmp(cfg.marker,'highlights')
  warning('using cfg.marker option -highlights- is no longer used, please use cfg.highlight')
  cfg.marker = 'off';
end

% check colormap is proper format and set it
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('topoplot(): Colormap must be a n x 3 matrix'); end
  colormap(cfg.colormap);
  ncolors = size(cfg.colormap,1);
else
  ncolors =[]; % let the low-level function deal with this
end

dtype  = ft_datatype(data);

% identify the interpretation of the functional data
switch dtype
  case 'raw'
    data   = ft_checkdata(data, 'datatype', 'timelock');
    dtype  = ft_datatype(data);
    dimord = data.dimord;
  case  {'timelock' 'freq' 'chan' 'unknown'}
    dimord = data.dimord;
  case 'comp'
    dimord = 'chan_comp';
  otherwise
end
dimtok = tokenize(dimord, '_');

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
    % Add a pseudo-axis with the component numbers:
    data.comp = 1:size(data.topo,2);
    % Specify the components
    if ~isempty(cfg.component)
      data.comp = cfg.component;
      data.topo = data.topo(:,cfg.component);
    end
    % Rename the field with topographic label information:
    data.label = data.topolabel;
    xparam = 'comp';
    yparam = '';
    cfg.parameter = ft_getopt(cfg, 'parameter', 'topo');
  otherwise
    % if the input data is not one of the standard data types, or if
    % the functional data is just one value per channel
    % in this case xparam, yparam are not defined
    % and the user should define the parameter
    if ~isfield(data, 'label'), error('the input data should at least contain a label-field'); end
    if ~isfield(cfg, 'parameter'), error('the configuration should at least contain a ''parameter'' field'); end
    if ~isfield(cfg, 'xparam'),
      cfg.xlim   = [1 1];
      xparam = '';
    end
end

if isfield(cfg, 'parameter') && ~isfield(data, cfg.parameter)
  error('cfg.parameter=%s is not present in data structure', cfg.parameter);
end

% user specified own fields, but no yparam (which is not asked in help)
if exist('xparam', 'var') && isfield(cfg, 'parameter') && ~exist('yparam', 'var')
  yparam = '';
end

% check whether rpt/subj is present and remove if necessary and whether
hasrpt = any(ismember(dimtok, {'rpt' 'subj'}));
if strcmp(dtype, 'timelock') && hasrpt,
  if ~isfield(data, cfg.parameter) || strcmp(cfg.parameter, 'individual')
    tmpcfg        = [];
    tmpcfg.trials = cfg.trials;
    data          = ft_timelockanalysis(tmpcfg, data);
    if ~strcmp(cfg.parameter, 'avg')
      % rename avg back into the parameter
      data.(cfg.parameter) = data.avg;
      data                 = rmfield(data, 'avg');
    end
    dimord        = data.dimord;
    dimtok        = tokenize(dimord, '_');
  else
    fprintf('input data contains repetitions, ignoring these and using ''%s'' field\n', cfg.parameter);
  end
elseif strcmp(dtype, 'freq') && hasrpt,
  % this also deals with fourier-spectra in the input
  % or with multiple subjects in a frequency domain stat-structure
  % on the fly computation of coherence spectrum is not supported
  if isfield(data, 'crsspctrm'), data = rmfield(data, 'crsspctrm'); end
  tmpcfg           = [];
  tmpcfg.trials    = cfg.trials;
  tmpcfg.jackknife = 'no';
  if isfield(cfg, 'parameter') && ~strcmp(cfg.parameter,'powspctrm')
    % freqdesctiptives will only work on the powspctrm field
    % hence a temporary copy of the data is needed
    tempdata.dimord    = data.dimord;
    tempdata.freq      = data.freq;
    tempdata.label     = data.label;
    tempdata.powspctrm = data.(cfg.parameter);
    if isfield(data, 'cfg'), tempdata.cfg = data.cfg; end
    tempdata           = ft_freqdescriptives(tmpcfg, tempdata);
    data.(cfg.parameter)  = tempdata.powspctrm;
    clear tempdata
  else
    data = ft_freqdescriptives(tmpcfg, data);
  end
  dimord = data.dimord;
  dimtok = tokenize(dimord, '_');
end

if isfield(data, 'label')
  cfg.channel = ft_channelselection(cfg.channel, data.label);
elseif isfield(data, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
end

% perform channel selection but only allow this when cfg.interactive = 'no'
if isfield(data, 'label') && strcmp(cfg.interactive, 'no')
  selchannel = ft_channelselection(cfg.channel, data.label);
elseif isfield(data, 'labelcmb') && strcmp(cfg.interactive, 'no')
  selchannel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
end

% Create time-series of small topoplots:
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

% Apply baseline correction:
if ~strcmp(cfg.baseline, 'no')
  if strcmp(xparam, 'freq') || strcmp(yparam, 'freq')
    data = ft_freqbaseline(cfg, data);
  else
    data = ft_timelockbaseline(cfg, data);
  end
end

% Handle the bivariate case

% Check for bivariate metric with 'chan_chan' in the dimord:
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;

% Check for bivariate metric with a labelcmb field:
haslabelcmb = isfield(data, 'labelcmb');

if (isfull || haslabelcmb) && (isfield(data, cfg.parameter) && ~strcmp(cfg.parameter, 'powspctrm'))
  % A reference channel is required:
  if ~isfield(cfg, 'refchannel')
    error('no reference channel is specified');
  end
  
  % check for refchannel being part of selection
  if ~strcmp(cfg.refchannel,'gui')
    if haslabelcmb
      cfg.refchannel = ft_channelselection(cfg.refchannel, unique(data.labelcmb(:)));
    else
      cfg.refchannel = ft_channelselection(cfg.refchannel, data.label);
    end
    if (isfull      && ~any(ismember(data.label, cfg.refchannel))) || ...
        (haslabelcmb && ~any(ismember(data.labelcmb(:), cfg.refchannel)))
      error('cfg.refchannel is a not present in the (selected) channels)')
    end
  end
  
  % Interactively select the reference channel
  if strcmp(cfg.refchannel, 'gui')
    % Open a single figure with the channel layout, the user can click on a reference channel
    h = clf;
    ft_plot_lay(cfg.layout, 'box', false);
    title('Select the reference channel by dragging a selection window, more than 1 channel can be selected...');
    % add the channel information to the figure
    info       = guidata(gcf);
    info.x     = cfg.layout.pos(:,1);
    info.y     = cfg.layout.pos(:,2);
    info.label = cfg.layout.label;
    guidata(h, info);
    % attach data to the figure with the current axis handle as a name
    dataname = fixname(num2str(double(gca)));
    setappdata(gcf,dataname,data);
    %set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'callback', {@select_topoplotER, cfg, data}});
    set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@select_topoplotER, cfg}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@select_topoplotER, cfg}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_topoplotER, cfg}, 'event', 'WindowButtonMotionFcn'});
    return
  end
  
  if ~isfull,
    % Convert 2-dimensional channel matrix to a single dimension:
    if isempty(cfg.directionality)
      sel1 = find(strcmp(cfg.refchannel, data.labelcmb(:,2)));
      sel2 = find(strcmp(cfg.refchannel, data.labelcmb(:,1)));
    elseif strcmp(cfg.directionality, 'outflow')
      sel1 = [];
      sel2 = find(strcmp(cfg.refchannel, data.labelcmb(:,1)));
    elseif strcmp(cfg.directionality, 'inflow')
      sel1 = find(strcmp(cfg.refchannel, data.labelcmb(:,2)));
      sel2 = [];
    end
    fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.parameter);
    if length(sel1)+length(sel2)==0
      error('there are no channels selected for plotting: you may need to look at the specification of cfg.directionality');
    end
    data.(cfg.parameter) = data.(cfg.parameter)([sel1;sel2],:,:);
    data.label     = [data.labelcmb(sel1,1);data.labelcmb(sel2,2)];
    data.labelcmb  = data.labelcmb([sel1;sel2],:);
    data           = rmfield(data, 'labelcmb');
  else
    % General case
    sel               = match_str(data.label, cfg.refchannel);
    siz               = [size(data.(cfg.parameter)) 1];
    if strcmp(cfg.directionality, 'inflow') || isempty(cfg.directionality)
      %the interpretation of 'inflow' and 'outflow' depend on
      %the definition in the bivariate representation of the data
      %in FieldTrip the row index 'causes' the column index channel
      %data.(cfg.parameter) = reshape(mean(data.(cfg.parameter)(:,sel,:),2),[siz(1) 1 siz(3:end)]);
      sel1 = 1:siz(1);
      sel2 = sel;
      meandir = 2;
    elseif strcmp(cfg.directionality, 'outflow')
      %data.(cfg.parameter) = reshape(mean(data.(cfg.parameter)(sel,:,:),1),[siz(1) 1 siz(3:end)]);
      sel1 = sel;
      sel2 = 1:siz(1);
      meandir = 1;
      
    elseif strcmp(cfg.directionality, 'inflow-outflow')
      % do the subtraction and recursively call the function again
      tmpcfg = cfg;
      tmpcfg.directionality = 'inflow';
      tmpdata = data;
      tmp     = data.(tmpcfg.parameter);
      siz     = [size(tmp) 1];
      for k = 1:siz(3)
        for m = 1:siz(4)
          tmp(:,:,k,m) = tmp(:,:,k,m)-tmp(:,:,k,m)';
        end
      end
      tmpdata.(tmpcfg.parameter) = tmp;
      ft_topoplotTFR(tmpcfg, tmpdata);
      return;
      
    elseif strcmp(cfg.directionality, 'outflow-inflow')
      % do the subtraction and recursively call the function again
      tmpcfg = cfg;
      tmpcfg.directionality = 'outflow';
      tmpdata = data;
      tmp     = data.(tmpcfg.parameter);
      siz     = [size(tmp) 1];
      for k = 1:siz(3)
        for m = 1:siz(4)
          tmp(:,:,k,m) = tmp(:,:,k,m)-tmp(:,:,k,m)';
        end
      end
      tmpdata.(tmpcfg.parameter) = tmp;
      ft_topoplotTFR(tmpcfg, tmpdata);
      return;
      
    end
  end
end

% Get physical min/max range of x:
if strcmp(cfg.xlim,'maxmin')
  xmin = min(data.(xparam));
  xmax = max(data.(xparam));
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Replace value with the index of the nearest bin
if ~isempty(xparam)
  xmin = nearest(data.(xparam), xmin);
  xmax = nearest(data.(xparam), xmax);
end

% Get physical min/max range of y:
if ~isempty(yparam)
  if strcmp(cfg.ylim,'maxmin')
    ymin = min(data.(yparam));
    ymax = max(data.(yparam));
  else
    ymin = cfg.ylim(1);
    ymax = cfg.ylim(2);
  end
  
  % Replace value with the index of the nearest bin:
  ymin = nearest(data.(yparam), ymin);
  ymax = nearest(data.(yparam), ymax);
end

% Take subselection of channels, this only works
% if the interactive mode is switched off
if exist('selchannel', 'var')
  sellab = match_str(data.label, selchannel);
  label  = data.label(sellab);
else
  sellab = 1:numel(data.label);
  label  = data.label;
end

if isfull
  sel1 = intersect(sel1, sellab);
  sel2 = intersect(sel2, sellab);
end

if ~isempty(cfg.parameter)
  % Make data vector with one value for each channel
  dat = data.(cfg.parameter);
  % get dimord dimensions
  ydim = find(strcmp(yparam, dimtok));
  xdim = find(strcmp(xparam, dimtok));
  zdim = setdiff(1:ndims(dat), [ydim xdim]);
  % and permute
  dat = permute(dat, [zdim(:)' ydim xdim]);
  
  if ~isempty(yparam)
    % time-frequency data
    if isfull
      dat = dat(sel1, sel2, ymin:ymax, xmin:xmax);
      dat = nanmean(nanmean(nanmean(dat, meandir), 4), 3);
    elseif haslabelcmb
      dat = dat(sellab, ymin:ymax, xmin:xmax);
      dat = nanmean(nanmean(dat, 3), 2);
    else
      dat = dat(sellab, ymin:ymax, xmin:xmax);
      dat = nanmean(nanmean(dat, 3), 2);
    end
  elseif ~isempty(cfg.component)
    % component data, nothing to do
  else
    % time or frequency data
    if isfull
      dat = dat(sel1, sel2, xmin:xmax);
      dat = nanmean(nanmean(dat, meandir), 3);
    elseif haslabelcmb
      dat = dat(sellab, xmin:xmax);
      dat = nanmean(dat, 2);
    else
      dat = dat(sellab, xmin:xmax);
      dat = nanmean(dat, 2);
    end
  end
  dat = dat(:);
  
else
  error('cannot make selection of data');
end

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
    if isfull
      msk = msk(sel1, sel2, ymin:ymax, xmin:xmax);
    elseif haslabelcmb
      msk = msk(sellab, ymin:ymax, xmin:xmax);
    else
      msk = msk(sellab, ymin:ymax, xmin:xmax);
    end
  elseif ~isempty(cfg.component)
    % component data, nothing to do
  else
    % time or frequency data
    if isfull
      msk = msk(sel1, sel2, xmin:xmax);
    elseif haslabelcmb
      msk = msk(sellab, xmin:xmax);
    else
      msk = msk(sellab, xmin:xmax);
    end
  end
  
  if size(msk,2)>1 || size(msk,3)>1
    warning('no masking possible for average over multiple latencies or frequencies -> cfg.maskparameter cleared')
    msk = [];
  end
  
else
  msk = [];
end


% Select the channels in the data that match with the layout:
[seldat, sellay] = match_str(label, cfg.layout.label);
if isempty(seldat)
  error('labels in data and labels in layout do not match');
end

dat = dat(seldat);
if ~isempty(msk)
  msk = msk(seldat);
end

% Select x and y coordinates and labels of the channels in the data
chanX = cfg.layout.pos(sellay,1);
chanY = cfg.layout.pos(sellay,2);
chanLabels = cfg.layout.label(sellay);

% Get physical min/max range of z:
if strcmp(cfg.zlim,'maxmin')
  zmin = min(dat);
  zmax = max(dat);
elseif strcmp(cfg.zlim,'maxabs')
  zmin = -max(max(abs(dat)));
  zmax = max(max(abs(dat)));
elseif strcmp(cfg.zlim,'zeromax')
  zmin = 0;
  zmax = max(dat);
elseif strcmp(cfg.zlim,'minzero')
  zmin = min(dat);
  zmax = 0;
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% make comment
if strcmp(cfg.comment, 'auto')
  comment = date;
  if ~isempty(xparam)
    if strcmp(cfg.xlim,'maxmin')
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, xparam, data.(xparam)(xmin), data.(xparam)(xmax));
    else
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, xparam, data.(xparam)(xmin), data.(xparam)(xmax));
    end
  end
  if ~isempty(yparam)
    if strcmp(cfg.ylim,'maxmin')
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, yparam, data.(yparam)(ymin), data.(yparam)(ymax));
    else
      comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, yparam, data.(yparam)(ymin), data.(yparam)(ymax));
    end
  end
  if ~isempty(cfg.parameter)
    comment = sprintf('%0s\n%0s=[%.3g %.3g]', comment, cfg.parameter, zmin, zmax);
  end
  if isfield(cfg,'refchannel')
    if iscell(cfg.refchannel)
      cfg.comment = sprintf('%s\nreference=%s %s', comment, cfg.refchannel{:});
    else
      cfg.comment = sprintf('%s\nreference=%s %s', comment, cfg.refchannel);
    end
  end
  cfg.comment = comment;
elseif strcmp(cfg.comment, 'xlim')
  if strcmp(cfg.xlim,'maxmin')
    comment = sprintf('%0s=[%.3g %.3g]', xparam, data.(xparam)(xmin), data.(xparam)(xmax));
  else
    comment = sprintf('%0s=[%.3g %.3g]', xparam, data.(xparam)(xmin), data.(xparam)(xmax));
  end
  if isfield(cfg,'refchannel')
    if iscell(cfg.refchannel)
      cfg.comment = sprintf('%s\nreference=%s %s', comment, cfg.refchannel{:});
    else
      cfg.comment = sprintf('%s\nreference=%s %s', comment, cfg.refchannel);
    end
  end
  cfg.comment = comment;
elseif ~ischar(cfg.comment)
  error('cfg.comment must be string');
end


% Specify the x and y coordinates of the comment
if strcmp(cfg.commentpos,'layout')
  ind_comment = strmatch('COMNT', cfg.layout.label);
  x_comment = cfg.layout.pos(ind_comment,1);
  y_comment = cfg.layout.pos(ind_comment,2);
elseif strcmp(cfg.commentpos,'lefttop')
  x_comment = -0.7;
  y_comment =  0.6;
  HorAlign = 'left';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos,'leftbottom')
  x_comment = -0.6;
  y_comment = -0.6;
  HorAlign = 'left';
  VerAlign = 'bottom';
elseif strcmp(cfg.commentpos,'middletop')
  x_comment =  0;
  y_comment =  0.75;
  HorAlign = 'center';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos,'middlebottom')
  x_comment =  0;
  y_comment = -0.7;
  HorAlign = 'center';
  VerAlign = 'bottom';
elseif strcmp(cfg.commentpos,'righttop')
  x_comment =  0.65;
  y_comment =  0.6;
  HorAlign = 'right';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos,'rightbottom')
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

% Set ft_plot_topo specific options
if strcmp(cfg.interplimits,'head')
  interplimits = 'mask';
else
  interplimits = cfg.interplimits;
end
if strcmp(cfg.style,'both');            style = 'surfiso';     end
if strcmp(cfg.style,'straight');        style = 'surf';        end
if strcmp(cfg.style,'contour');         style = 'iso';         end
if strcmp(cfg.style,'fill');            style = 'isofill';     end
if strcmp(cfg.style,'straight_imsat');  style = 'imsat';       end
if strcmp(cfg.style,'both_imsat');      style = 'imsatiso';    end

% check for nans
nanInds = isnan(dat);
if strcmp(cfg.interpolatenan,'yes') && any(nanInds)
  warning('removing NaNs from the data');
  chanX(nanInds) = [];
  chanY(nanInds) = [];
  dat(nanInds)   = [];
  if ~isempty(msk)
    msk(nanInds) = [];
  end
end

% Draw plot
if ~strcmp(cfg.style,'blank')
  opt = {'interpmethod',cfg.interpolation,...
    'interplim',interplimits,...
    'gridscale',cfg.gridscale,...
    'outline',cfg.layout.outline,...
    'shading',cfg.shading,...
    'isolines',cfg.contournum,...
    'mask',cfg.layout.mask,...
    'style',style,...
    'datmask', msk};
  if strcmp(style,'imsat') || strcmp(style,'imsatiso')
    % add clim to opt
    opt = [opt {'clim',[zmin zmax],'ncolors',ncolors}];
  end
  ft_plot_topo(chanX,chanY,dat,opt{:});
elseif ~strcmp(cfg.style,'blank')
  ft_plot_lay(cfg.layout,'box','no','label','no','point','no')
end

% Plotting markers for channels and/or highlighting a selection of channels
highlightchansel = []; % used for remembering selection of channels
templay.outline = cfg.layout.outline;
templay.mask    = cfg.layout.mask;
% For Highlight (channel-selection)
for icell = 1:length(cfg.highlight)
  if ~strcmp(cfg.highlight{icell},'off')
    [dum labelindex]   = match_str(ft_channelselection(cfg.highlightchannel{icell}, data.label), cfg.layout.label);
    highlightchansel   = [highlightchansel; match_str(data.label,ft_channelselection(cfg.highlightchannel{icell}, data.label))];
    templay.pos        = cfg.layout.pos(labelindex,:);
    templay.width      = cfg.layout.width(labelindex);
    templay.height     = cfg.layout.height(labelindex);
    templay.label      = cfg.layout.label(labelindex);
    if strcmp(cfg.highlight{icell}, 'labels') || strcmp(cfg.highlight{icell}, 'numbers')
      labelflg = 1;
    else
      labelflg = 0;
    end
    if strcmp(cfg.highlight{icell}, 'numbers')
      for ichan = 1:length(labelindex)
        templay.label{ichan} = num2str(match_str(data.label,templay.label{ichan}));
      end
    end
    ft_plot_lay(templay,'box','no','label',labelflg,'point','yes',...
      'pointsymbol',cfg.highlightsymbol{icell},...
      'pointcolor',cfg.highlightcolor{icell},...
      'pointsize',cfg.highlightsize{icell},...
      'fontsize',cfg.highlightfontsize{icell},...
      'labeloffset',cfg.labeloffset)
  end
end % for icell
% For Markers (all channels)
if ~strcmp(cfg.marker,'off')
  channelsToMark = 1:length(data.label);
  if strcmp(cfg.interpolatenan,'no')
    channelsNotMark = highlightchansel;
  else
    channelsNotMark = union(find(isnan(dat)),highlightchansel);
  end
  channelsToMark(channelsNotMark) = [];
  [dum labelindex] = match_str(ft_channelselection(channelsToMark, data.label),cfg.layout.label);
  templay.pos      = cfg.layout.pos(labelindex,:);
  templay.width    = cfg.layout.width(labelindex);
  templay.height   = cfg.layout.height(labelindex);
  templay.label    = cfg.layout.label(labelindex);
  if strcmp(cfg.marker, 'labels') || strcmp(cfg.marker, 'numbers')
    labelflg = 1;
  else
    labelflg = 0;
  end
  if strcmp(cfg.marker, 'numbers')
    for ichan = 1:length(labelindex)
      templay.label{ichan} = num2str(match_str(data.label,templay.label{ichan}));
    end
  end
  ft_plot_lay(templay,'box','no','label',labelflg,'point','yes',...
    'pointsymbol',cfg.markersymbol,...
    'pointcolor',cfg.markercolor,...
    'pointsize',cfg.markersize,...
    'fontsize',cfg.markerfontsize,...
    'labeloffset',cfg.labeloffset)
end


if(isfield(cfg,'vector'))
    vecX = mean(real(data.(cfg.vector)(:,xmin:xmax)),2);
    vecY = mean(imag(data.(cfg.vector)(:,xmin:xmax)),2);
    
    % scale quiver relative to largest gradiometer sample
    k=0.15/max([max(abs(real(data.(cfg.vector)(:)))) max(abs(imag(data.(cfg.vector)(:))))]);
    quiver(chanX, chanY, k*vecX, k*vecY,0,'red');
end

% Write comment
if ~strcmp(cfg.comment,'no')
  if strcmp(cfg.commentpos, 'title')
    title(cfg.comment, 'Fontsize', cfg.fontsize);
  else
    ft_plot_text(x_comment,y_comment, cfg.comment, 'Fontsize', cfg.fontsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Fontweight', cfg.fontweight);
  end
end

% plot colorbar:
if isfield(cfg, 'colorbar')
  if strcmp(cfg.colorbar, 'yes')
    colorbar;
  elseif ~strcmp(cfg.colorbar, 'no')
    colorbar('location',cfg.colorbar);
  end
end

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

% The remainder of the code is meant to make the figure interactive
hold on;

% Set colour axis
if ~strcmp(cfg.style, 'blank')
  caxis([zmin zmax]);
end

if strcmp('yes',cfg.hotkeys)
  %  Attach data and cfg to figure and attach a key listener to the figure
  set(gcf, 'KeyPressFcn', {@key_sub, zmin, zmax})
end

% Make the figure interactive
if strcmp(cfg.interactive, 'yes')
    % add the cfg/data/channel information to the figure under identifier linked to this axis
    ident                    = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
    set(gca,'tag',ident);
    info                     = guidata(gcf);
    info.(ident).x           = cfg.layout.pos(:, 1);
    info.(ident).y           = cfg.layout.pos(:, 2);
    info.(ident).label       = cfg.layout.label;
    info.(ident).cfg         = cfg;
    info.(ident).datvarargin = varargin(1:Ndata);
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
    warning('unsupported dimord "%s" for interactive plotting', data.dimord);
  end
end

% set the figure window title, but only if the user has not changed it
if isempty(get(gcf, 'Name'))
  if isfield(cfg,'funcname')
    funcname = cfg.funcname;
  else
    funcname = mfilename;
  end
  if isfield(cfg,'dataname')
    if iscell(cfg.dataname)
      dataname = cfg.dataname{indx};
    else
      dataname = cfg.dataname;
    end
  elseif nargin > 1
    dataname = {inputname(2)};
    for k = 2:Ndata
      dataname{end+1} = inputname(k+1);
    end
  else % data provided through cfg.inputfile
    dataname = cfg.inputfile;
  end
  
  if isempty(cfg.figurename)
    dataname_str=join_str(', ', dataname);

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
% SUBFUNCTION which is called after selecting channels in case of cfg.refchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotER(label, cfg, varargin)
% get appdata belonging to current axis
dataname = fixname(num2str(double(gca)));
data = getappdata(gcf, dataname);

if isfield(cfg, 'inputfile')
  % the reading has already been done and varargin contains the data
  cfg = rmfield(cfg, 'inputfile');
end
cfg.refchannel = label;
fprintf('selected cfg.refchannel = ''%s''\n', cfg.refchannel{:});
cfg.highlight = 'on';
cfg.highlightsymbol  = '.';
cfg.highlightcolor   = 'r';
cfg.highlightsize = 20;
cfg.highlightchannel =  cfg.refchannel;
% ensure that the new figure appears at the same position
f = figure('Position', get(gcf, 'Position'), 'Visible', get(gcf, 'Visible'));
ft_topoplotER(cfg, data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label)
% fetch cfg/data based on axis indentifier given as tag
ident       = get(gca,'tag');
info        = guidata(gcf);
cfg         = info.(ident).cfg;
datvarargin = info.(ident).datvarargin;
if ~isempty(label)
  
  if isfield(cfg, 'inputfile')
    % the reading has already been done and varargin contains the data
    cfg = rmfield(cfg, 'inputfile');
  end
  cfg.xlim = 'maxmin';
  cfg.channel = label;
  % if user specified a zlim, copy it over to the ylim of singleplot
  if isfield(cfg, 'zlim')
    cfg.ylim = cfg.zlim;
    cfg = rmfield(cfg, 'zlim');
  end
  fprintf('selected cfg.channel = {');
  for i=1:(length(cfg.channel)-1)
    fprintf('''%s'', ', cfg.channel{i});
  end
  fprintf('''%s''}\n', cfg.channel{end});
  % ensure that the new figure appears at the same position
  f = figure('Position', get(gcf, 'Position'), 'Visible', get(gcf, 'Visible'));
  ft_singleplotER(cfg, datvarargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotTFR(label)
ident       = get(gca,'tag');
info        = guidata(gcf);
cfg         = info.(ident).cfg;
datvarargin = info.(ident).datvarargin;
if ~isempty(label)
  
  if isfield(cfg, 'inputfile')
    % the reading has already been done and varargin contains the data
    cfg = rmfield(cfg, 'inputfile');
  end
  cfg.xlim = 'maxmin';
  cfg.ylim = 'maxmin';
  cfg.channel = label;
  fprintf('selected cfg.channel = {');
  for i=1:(length(cfg.channel)-1)
    fprintf('''%s'', ', cfg.channel{i});
  end
  fprintf('''%s''}\n', cfg.channel{end});
  % ensure that the new figure appears at the same position
  f = figure('Position', get(gcf, 'Position'), 'Visible', get(gcf, 'Visible'));
  ft_singleplotTFR(cfg, datvarargin{:});
end


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
