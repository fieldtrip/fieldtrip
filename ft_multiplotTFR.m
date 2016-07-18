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
%   cfg.maskparameter    = field in the data to be used for opacity masking of data
%   cfg.maskstyle        = style used to masking, 'opacity', 'saturation' or 'outline' (default = 'opacity')
%                          use 'saturation' or 'outline' when saving to vector-format (like *.eps) to avoid all
%                          sorts of image-problems (currently only possible with a white backgroud)
%   cfg.maskalpha        = alpha value between 0 (transparant) and 1 (opaque) used for masking areas dictated by cfg.maskparameter (default = 1)
%   cfg.masknans         = 'yes' or 'no' (default = 'yes')
%   cfg.xlim             = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim             = 'maxmin' or [ymin ymax] (default = 'maxmin')
%   cfg.zlim             = plotting limits for color dimension, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.gradscale        = number, scaling to apply to the MEG gradiometer channels prior to display
%   cfg.magscale         = number, scaling to apply to the MEG magnetometer channels prior to display
%   cfg.channel          = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.refchannel       = name of reference channel for visualising connectivity, can be 'gui'
%   cfg.baseline         = 'yes', 'no' or [time1 time2] (default = 'no'), see FT_FREQBASELINE
%   cfg.baselinetype     = 'absolute', 'relative', 'relchange' or 'db' (default = 'absolute')
%   cfg.trials           = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.box              = 'yes', 'no' (default = 'no' if maskparameter given default = 'yes')
%                          Draw a box around each graph
%   cfg.hotkeys          = enables hotkeys (up/down arrows) for dynamic colorbar adjustment
%   cfg.colorbar         = 'yes', 'no' (default = 'no')
%   cfg.colormap         = any sized colormap, see COLORMAP
%   cfg.comment          = string of text (default = date + zlimits)
%                          Add 'comment' to graph (according to COMNT in the layout)
%   cfg.showlabels       = 'yes', 'no' (default = 'no')
%   cfg.showoutline      = 'yes', 'no' (default = 'no')
%   cfg.fontsize         = font size of comment and labels (if present) (default = 8)
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
%                         the supported ways (see below).
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
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure. For this particular function, the
% data should be provided as a cell array.
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
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'freq');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'unused', {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamed', {'matrixside', 'directionality'});
cfg = ft_checkconfig(cfg, 'renamed', {'cohrefchannel', 'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed', {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim', 'absmax', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedback', 'inflow'});
cfg = ft_checkconfig(cfg, 'deprecated', {'xparam', 'yparam'});

% set the defaults
cfg.baseline       = ft_getopt(cfg, 'baseline', 'no');
cfg.baselinetype   = ft_getopt(cfg, 'baselinetype', 'absolute');
cfg.trials         = ft_getopt(cfg, 'trials', 'all', 1);
cfg.xlim           = ft_getopt(cfg, 'xlim', 'maxmin');
cfg.ylim           = ft_getopt(cfg, 'ylim', 'maxmin');
cfg.zlim           = ft_getopt(cfg, 'zlim', 'maxmin');
cfg.magscale       = ft_getopt(cfg, 'magscale', 1);
cfg.gradscale      = ft_getopt(cfg, 'gradscale', 1);
cfg.colorbar       = ft_getopt(cfg, 'colorbar', 'no');
cfg.comment        = ft_getopt(cfg, 'comment', date);
cfg.showlabels     = ft_getopt(cfg, 'showlabels', 'no');
cfg.showoutline    = ft_getopt(cfg, 'showoutline', 'no');
cfg.channel        = ft_getopt(cfg, 'channel', 'all');
cfg.fontsize       = ft_getopt(cfg, 'fontsize', 8);
cfg.interactive    = ft_getopt(cfg, 'interactive', 'yes');
cfg.hotkeys        = ft_getopt(cfg, 'hotkeys', 'no');
cfg.renderer       = ft_getopt(cfg, 'renderer'); % let MATLAB decide on default
cfg.orient         = ft_getopt(cfg, 'orient', 'landscape');
cfg.maskalpha      = ft_getopt(cfg, 'maskalpha', 1);
cfg.masknans       = ft_getopt(cfg, 'masknans', 'yes');
cfg.maskparameter  = ft_getopt(cfg, 'maskparameter');
cfg.maskstyle      = ft_getopt(cfg, 'maskstyle', 'opacity');
cfg.directionality = ft_getopt(cfg, 'directionality', '');
cfg.figurename     = ft_getopt(cfg, 'figurename');
if ~isfield(cfg, 'box')
  if ~isempty(cfg.maskparameter)
    cfg.box = 'yes';
  else
    cfg.box = 'no';
  end
end
if numel(findobj(gcf, 'type', 'axes', '-not', 'tag', 'ft-colorbar')) > 1 && strcmp(cfg.interactive, 'yes')
  warning('using cfg.interactive = ''yes'' in subplots is not supported, setting cfg.interactive = ''no''')
  cfg.interactive = 'no';
end

dimord = data.dimord;
dimtok = tokenize(dimord, '_');

% Set x/y/parameter defaults
if ~any(ismember(dimtok, 'time'))
  error('input data needs a time dimension');
else
  xparam = 'time';
  yparam = 'freq';
  cfg.parameter = ft_getopt(cfg, 'parameter', 'powspctrm');
end

if isfield(cfg, 'channel') && isfield(data, 'label')
  cfg.channel = ft_channelselection(cfg.channel, data.label);
elseif isfield(cfg, 'channel') && isfield(data, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
end

% perform channel selection but only allow this when cfg.interactive = 'no'
if isfield(data, 'label') && strcmp(cfg.interactive, 'no')
  selchannel = ft_channelselection(cfg.channel, data.label);
elseif isfield(data, 'labelcmb') && strcmp(cfg.interactive, 'no')
  selchannel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
end

% check whether rpt/subj is present and remove if necessary and whether
hasrpt = any(ismember(dimtok, {'rpt' 'subj'}));
if hasrpt,
  % this also deals with fourier-spectra in the input
  % or with multiple subjects in a frequency domain stat-structure
  % on the fly computation of coherence spectrum is not supported
  if isfield(data, 'crsspctrm'),
    data = rmfield(data, 'crsspctrm');
  end
  % keep mask-parameter if it is set
  if ~isempty(cfg.maskparameter)
    tempmask = data.(cfg.maskparameter);
  end
  tmpcfg           = [];
  tmpcfg.trials    = cfg.trials;
  tmpcfg.jackknife = 'no';
  if isfield(cfg, 'parameter') && ~strcmp(cfg.parameter, 'powspctrm')
    % freqdesctiptives will only work on the powspctrm field
    % hence a temporary copy of the data is needed
    tempdata.dimord    = data.dimord;
    tempdata.freq      = data.freq;
    tempdata.label     = data.label;
    tempdata.powspctrm = data.(cfg.parameter);
    if isfield(data, 'cfg') tempdata.cfg = data.cfg; end
    tempdata           = ft_freqdescriptives(tmpcfg, tempdata);
    data.(cfg.parameter)  = tempdata.powspctrm;
    clear tempdata
  else
    data = ft_freqdescriptives(tmpcfg, data);
  end
  % put mask-parameter back if it is set
  if ~isempty(cfg.maskparameter)
    data.(cfg.maskparameter) = tempmask;
  end
  dimord = data.dimord;
  dimtok = tokenize(dimord, '_');
end % if hasrpt

% Read or create the layout that will be used for plotting:
cla;
hold on
lay = ft_prepare_layout(cfg, data);
cfg.layout = lay;


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

% Handle the bivariate case

% Check for bivariate metric with 'chan_chan' in the dimord
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;

% Check for bivariate metric with a labelcmb
haslabelcmb = isfield(data, 'labelcmb');

if (isfull || haslabelcmb) && (isfield(data, cfg.parameter) && ~strcmp(cfg.parameter, 'powspctrm'))
  % A reference channel is required:
  if ~isfield(cfg, 'refchannel')
    error('no reference channel is specified');
  end

  % check for refchannel being part of selection
  if ~strcmp(cfg.refchannel, 'gui')
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
    ft_plot_lay(lay, 'box', false);
    title('Select the reference channel by dragging a selection window, more than 1 channel can be selected...');
    % add the channel information to the figure
    info       = guidata(gcf);
    info.x     = lay.pos(:, 1);
    info.y     = lay.pos(:, 2);
    info.label = lay.label;
    info.dataname = '';
    guidata(h, info);
    %set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'callback', {@select_topoplotER, cfg, data}});
    set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
    return
  end

  if ~isfull,
    % Convert 2-dimensional channel matrix to a single dimension:
    if isempty(cfg.directionality)
      sel1 = find(strcmp(cfg.refchannel, data.labelcmb(:, 2)));
      sel2 = find(strcmp(cfg.refchannel, data.labelcmb(:, 1)));
    elseif strcmp(cfg.directionality, 'outflow')
      sel1 = [];
      sel2 = find(strcmp(cfg.refchannel, data.labelcmb(:, 1)));
    elseif strcmp(cfg.directionality, 'inflow')
      sel1 = find(strcmp(cfg.refchannel, data.labelcmb(:, 2)));
      sel2 = [];
    end
    fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.parameter);
    if length(sel1)+length(sel2)==0
      error('there are no channels selected for plotting: you may need to look at the specification of cfg.directionality');
    end
    data.(cfg.parameter) = data.(cfg.parameter)([sel1;sel2], :, :);
    data.label     = [data.labelcmb(sel1, 1);data.labelcmb(sel2, 2)];
    data.labelcmb  = data.labelcmb([sel1;sel2], :);
    %data           = rmfield(data, 'labelcmb');
  else
    % General case
    sel               = match_str(data.label, cfg.refchannel);
    siz               = [size(data.(cfg.parameter)) 1];
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
      error('cfg.directionality = ''ff-fd'' is not supported anymore, you have to manually subtract the two before the call to ft_multiplotTFR');
    elseif strcmp(cfg.directionality, 'fd-ff')
      error('cfg.directionality = ''fd-ff'' is not supported anymore, you have to manually subtract the two before the call to ft_multiplotTFR');
    end %if directionality
  end %if ~isfull
end %handle the bivariate data


% Get physical x-axis range:
if strcmp(cfg.xlim, 'maxmin')
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

% Get physical y-axis range:
if strcmp(cfg.ylim, 'maxmin')
  ymin = min(data.(yparam));
  ymax = max(data.(yparam));
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Replace value with the index of the nearest bin
if ~isempty(yparam)
  ymin = nearest(data.(yparam), ymin);
  ymax = nearest(data.(yparam), ymax);
end

% test if X and Y are linearly spaced (to within 10^-12): % FROM UIMAGE
x = data.(xparam)(xmin:xmax);
y = data.(yparam)(ymin:ymax);
dx = min(diff(x));  % smallest interval for X
dy = min(diff(y));  % smallest interval for Y
evenx = all(abs(diff(x)/dx-1)<1e-12);     % true if X is linearly spaced
eveny = all(abs(diff(y)/dy-1)<1e-12);     % true if Y is linearly spaced

if ~evenx || ~eveny
  warning('(one of the) axis is/are not evenly spaced, but plots are made as if axis are linear')
end

% Take subselection of channels, this only works
% in the interactive mode
if exist('selchannel', 'var')
  sellab = match_str(data.label, selchannel);
  label  = data.label(sellab);
else
  sellab = 1:numel(data.label);
  label  = data.label;
end

dat = data.(cfg.parameter);
% get dimord dimensions
dims = textscan(data.dimord, '%s', 'Delimiter', '_');
dims = dims{1};
ydim = find(strcmp(yparam, dims));
xdim = find(strcmp(xparam, dims));
zdim = setdiff(1:ndims(dat), [ydim xdim]);
% and permute
dat = permute(dat, [zdim(:)' ydim xdim]);
if isfull
  dat = dat(sel1, sel2, ymin:ymax, xmin:xmax);
  dat = nanmean(dat, meandir);
  siz = size(dat);
  dat = reshape(dat, [max(siz(1:2)) siz(3) siz(4)]);
  dat = dat(sellab, :, :);
% this makes no sense, so COMMENTED OUT AS OF FEBURARY 22 2012
% elseif haslabelcmb
%   dat = dat(sellab, ymin:ymax, xmin:xmax);
else
  dat = dat(sellab, ymin:ymax, xmin:xmax);
end

if ~isempty(cfg.maskparameter)
  mask = data.(cfg.maskparameter);
  mask = permute(mask, [zdim(:)' ydim xdim]);
  if isfull && cfg.maskalpha == 1
    mask = mask(sel1, sel2, ymin:ymax, xmin:xmax);
    mask = nanmean(nanmean(nanmean(mask, meandir), 4), 3);
  elseif haslabelcmb && cfg.maskalpha == 1
    mask = mask(sellab, ymin:ymax, xmin:xmax);
    %mask = nanmean(nanmean(mask, 3), 2);
  elseif cfg.maskalpha == 1
    mask = mask(sellab, ymin:ymax, xmin:xmax);
    %mask = nanmean(nanmean(mask, 3), 2);
  elseif isfull && cfg.maskalpha ~= 1
    maskl = mask(sel1, sel2, ymin:ymax, xmin:xmax); %% check this for full representation
    mask = zeros(size(maskl));
    mask(maskl) = 1;
    mask(~maskl) = cfg.maskalpha;
  elseif haslabelcmb && cfg.maskalpha ~= 1
    maskl = mask(sellab, ymin:ymax, xmin:xmax);
    mask = zeros(size(maskl));
    mask(maskl) = 1;
    mask(~maskl) = cfg.maskalpha;
  elseif cfg.maskalpha ~= 1
    maskl = mask(sellab, ymin:ymax, xmin:xmax);
    mask = zeros(size(maskl));
    mask(maskl) = 1;
    mask(~maskl) = cfg.maskalpha;
  end
end

% Select the channels in the data that match with the layout:
[chanseldat, chansellay] = match_str(label, lay.label);
if isempty(chanseldat)
  error('labels in data and labels in layout do not match');
end

% if magnetometer/gradiometer scaling is requested, get indices for
% channels
if (cfg.magscale ~= 1)
  magInd = match_str(label, ft_channelselection('MEGMAG', label));
end
if (cfg.gradscale ~= 1)
  gradInd = match_str(label, ft_channelselection('MEGGRAD', label));
end

datsel = dat(chanseldat, :, :);
if ~isempty(cfg.maskparameter)
  maskdat = mask(chanseldat, :, :);
end

% Select x and y coordinates and labels of the channels in the data
chanX = lay.pos(chansellay, 1);
chanY = lay.pos(chansellay, 2);
chanWidth  = lay.width(chansellay);
chanHeight = lay.height(chansellay);

% Get physical z-axis range (color axis):
if strcmp(cfg.zlim, 'maxmin')
  zmin = min(datsel(:));
  zmax = max(datsel(:));
elseif strcmp(cfg.zlim, 'maxabs')
  zmin = -max(abs(datsel(:)));
  zmax = max(abs(datsel(:)));
elseif strcmp(cfg.zlim, 'zeromax')
  zmin = 0;
  zmax = max(datsel(:));
elseif strcmp(cfg.zlim, 'minzero')
  zmin = min(datsel(:));
  zmax = 0;
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% set colormap
if isfield(cfg, 'colormap')
  if size(cfg.colormap, 2)~=3, error('multiplotTFR(): Colormap must be a n x 3 matrix'); end
  set(gcf, 'colormap', cfg.colormap);
end

% Plot channels:
for k=1:length(chanseldat)
  % Get cdata:
  cdata = shiftdim(datsel(k, :, :));
  if ~isempty(cfg.maskparameter)
    mdata = shiftdim(maskdat(k, :, :));
  end

  % scale if needed
  if (cfg.magscale ~= 1 && any(magInd == chanseldat(k)))
    cdata = cdata .* cfg.magscale;
  end
  if (cfg.gradscale ~= 1 && any(gradInd == chanseldat(k)))
    cdata = cdata .* cfg.gradscale;
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
end % for chanseldat

% write comment:
k = cellstrmatch('COMNT', lay.label);
if ~isempty(k)
  comment = cfg.comment;
  comment = sprintf('%0s\nxlim=[%.3g %.3g]', comment, data.(xparam)(xmin), data.(xparam)(xmax));
  comment = sprintf('%0s\nylim=[%.3g %.3g]', comment, data.(yparam)(ymin), data.(yparam)(ymax));
  comment = sprintf('%0s\nzlim=[%.3g %.3g]', comment, zmin, zmax);
  ft_plot_text(lay.pos(k, 1), lay.pos(k, 2), sprintf(comment), 'Fontsize', cfg.fontsize);
end

% plot scale:
k = cellstrmatch('SCALE', lay.label);
if ~isempty(k)
  % Get average cdata across channels:
  cdata = shiftdim(mean(datsel, 1));

  % Draw plot (and mask Nan's with maskfield if requested)
  if isequal(cfg.masknans, 'yes') && isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = double(mask);
    ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask, 'hpos', lay.pos(k, 1), 'vpos', lay.pos(k, 2), 'width', lay.width(k, 1), 'height', lay.height(k, 1))
  elseif isequal(cfg.masknans, 'yes') && ~isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = mask .* mdata;
    mask = double(mask);
    ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask, 'hpos', lay.pos(k, 1), 'vpos', lay.pos(k, 2), 'width', lay.width(k, 1), 'height', lay.height(k, 1))
  elseif isequal(cfg.masknans, 'no') && ~isempty(cfg.maskparameter)
    mask = mdata;
    mask = double(mask);
    ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'highlightstyle', cfg.maskstyle, 'highlight', mask, 'hpos', lay.pos(k, 1), 'vpos', lay.pos(k, 2), 'width', lay.width(k, 1), 'height', lay.height(k, 1))
  else
    ft_plot_matrix(cdata, 'clim', [zmin zmax], 'tag', 'cip', 'hpos', lay.pos(k, 1), 'vpos', lay.pos(k, 2), 'width', lay.width(k, 1), 'height', lay.height(k, 1))
  end
  % Currently the handle isn't being used below, this is here for possible use in the future
  h = findobj('tag', 'cip');

end

% plot layout
boxflg     = istrue(cfg.box);
labelflg   = istrue(cfg.showlabels);
outlineflg = istrue(cfg.showoutline);
ft_plot_lay(lay, 'box', boxflg, 'label', labelflg, 'outline', outlineflg, 'point', 'no', 'mask', 'no');

% plot colorbar:
if isfield(cfg, 'colorbar') && (strcmp(cfg.colorbar, 'yes'))
  colorbar;
end

% Set colour axis
caxis([zmin zmax]);
if strcmp('yes', cfg.hotkeys)
  %  Attach data and cfg to figure and attach a key listener to the figure
  set(gcf, 'KeyPressFcn', {@key_sub, zmin, zmax})
end

% set the figure window title
if isempty(get(gcf, 'Name'))
  if isfield(cfg, 'funcname')
    funcname = cfg.funcname;
  else
    funcname = mfilename;
  end

  if isfield(cfg, 'dataname')
      dataname = cfg.dataname;
  elseif nargin > 1
    dataname = inputname(2);
  else % data provided through cfg.inputfile
    dataname = cfg.inputfile;
  end

  if isempty(cfg.figurename)
    set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), funcname, dataname));
    set(gcf, 'NumberTitle', 'off');
  else
    set(gcf, 'name', cfg.figurename);
    set(gcf, 'NumberTitle', 'off');
  end
else
  funcname = '';
  dataname = '';
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
    % add the channel information to the figure
    info       = guidata(gcf);
    info.x     = lay.pos(:, 1);
    info.y     = lay.pos(:, 2);
    info.label = lay.label;
    info.dataname = dataname;
    guidata(gcf, info);

    set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
end

axis tight
axis off
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
ft_postamble previous data
ft_postamble provenance

% add a menu to the figure
% also, delete any possibly existing previous menu, this is safe because delete([]) does nothing
delete(findobj(gcf, 'type', 'uimenu', 'label', 'FieldTrip'));
ftmenu = uimenu(gcf, 'Label', 'FieldTrip');
uimenu(ftmenu, 'Label', 'Show pipeline', 'Callback', {@menu_pipeline, cfg});
uimenu(ftmenu, 'Label', 'About', 'Callback', @menu_about);

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
% SUBFUNCTION which is called by ft_select_channel in case cfg.refchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_multiplotTFR(label, cfg, varargin)
if isfield(cfg, 'inputfile')
  % the reading has already been done and varargin contains the data
  cfg = rmfield(cfg, 'inputfile');
end

% put data name in here, this cannot be resolved by other means
info = guidata(gcf);
cfg.dataname = info.dataname;

cfg.refchannel = label;
fprintf('selected cfg.refchannel = ''%s''\n', join_str(', ', cfg.refchannel));
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_multiplotTFR(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotTFR(label, cfg, varargin)
if ~isempty(label)
  if isfield(cfg, 'inputfile')
    % the reading has already been done and varargin contains the data
    cfg = rmfield(cfg, 'inputfile');
  end
  cfg.channel = label;

  % make sure ft_singleplotTFR does not apply a baseline correction again
  cfg.baseline = 'no';

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
  ft_singleplotTFR(cfg, varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata, varargin)
incr = (max(caxis)-min(caxis)) /10;
% symmetrically scale color bar down by 10 percent
if strcmp(eventdata.Key, 'uparrow')
  caxis([min(caxis)-incr max(caxis)+incr]);
% symmetrically scale color bar up by 10 percent
elseif strcmp(eventdata.Key, 'downarrow')
  caxis([min(caxis)+incr max(caxis)-incr]);
% resort to minmax of data for colorbar
elseif strcmp(eventdata.Key, 'm')
  caxis([varargin{1} varargin{2}]);
end
