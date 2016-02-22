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
%   cfg.maskstyle      = style used to masking, 'opacity', 'saturation' or 'outline' (default = 'opacity')
%                        use 'saturation' or 'outline' when saving to vector-format (like *.eps) to avoid all sorts of image-problems
%   cfg.maskalpha      = alpha value between 0 (transparant) and 1 (opaque) used for masking areas dictated by cfg.maskparameter (default = 1)
%   cfg.masknans       = 'yes' or 'no' (default = 'yes')
%   cfg.xlim           = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim           = 'maxmin' or [ymin ymax] (default = 'maxmin')
%   cfg.zlim           = plotting limits for color dimension, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.baseline       = 'yes','no' or [time1 time2] (default = 'no'), see FT_FREQBASELINE
%   cfg.baselinetype   = 'absolute', 'relative', 'relchange' or 'db' (default = 'absolute')
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                        see FT_CHANNELSELECTION for details
%   cfg.refchannel     = name of reference channel for visualising connectivity, can be 'gui'
%   cfg.fontsize       = font size of title (default = 8)
%   cfg.hotkeys           = enables hotkeys (up/down arrows) for dynamic colorbar adjustment
%   cfg.colormap       = any sized colormap, see COLORMAP
%   cfg.colorbar       = 'yes', 'no' (default = 'yes')
%   cfg.interactive    = Interactive plot 'yes' or 'no' (default = 'yes')
%                        In a interactive plot you can select areas and produce a new
%                        interactive plot when a selected area is clicked. Multiple areas
%                        can be selected by holding down the SHIFT key.
%   cfg.renderer       = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
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

% Copyright (C) 2005-2006, F.C. Donders Centre
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

% Set the defaults:
cfg.baseline       = ft_getopt(cfg, 'baseline',     'no');
cfg.baselinetype   = ft_getopt(cfg, 'baselinetype', 'absolute');
cfg.trials         = ft_getopt(cfg, 'trials',       'all', 1);
cfg.xlim           = ft_getopt(cfg, 'xlim',         'maxmin');
cfg.ylim           = ft_getopt(cfg, 'ylim',         'maxmin');
cfg.zlim           = ft_getopt(cfg, 'zlim',         'maxmin');
cfg.fontsize       = ft_getopt(cfg, 'fontsize',      8);
cfg.colorbar       = ft_getopt(cfg, 'colorbar',     'yes');
cfg.interactive    = ft_getopt(cfg, 'interactive',  'yes');
cfg.hotkeys        = ft_getopt(cfg, 'hotkeys',      'no');
cfg.renderer       = ft_getopt(cfg, 'renderer',      []);
cfg.maskalpha      = ft_getopt(cfg, 'maskalpha',     1);
cfg.maskparameter  = ft_getopt(cfg, 'maskparameter', []);
cfg.maskstyle      = ft_getopt(cfg, 'maskstyle',    'opacity');
cfg.channel        = ft_getopt(cfg, 'channel',      'all');
cfg.masknans       = ft_getopt(cfg, 'masknans',     'yes');
cfg.directionality = ft_getopt(cfg, 'directionality',[]);
cfg.figurename     = ft_getopt(cfg, 'figurename',    []);
cfg.parameter      = ft_getopt(cfg, 'parameter', 'powspctrm');

dimord = getdimord(data, cfg.parameter);
dimtok = tokenize(dimord, '_');

% Set x/y/parameter defaults
if ~any(ismember(dimtok, 'time'))
  error('input data needs a time dimension');
else
  xparam = 'time';
  yparam = 'freq';
end

if isfield(cfg, 'channel') && isfield(data, 'label')
  cfg.channel = ft_channelselection(cfg.channel, data.label);
elseif isfield(cfg, 'channel') && isfield(data, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
end

if isempty(cfg.channel)
  error('no channels selected');
end

if ~isfield(data, cfg.parameter)
  error('data has no field ''%s''', cfg.parameter);
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

  tmpcfg           = [];
  tmpcfg.trials    = cfg.trials;
  tmpcfg.jackknife = 'no';
  % keep mask-parameter if it is set
  if ~isempty(cfg.maskparameter)
    tempmask = data.(cfg.maskparameter);
  end
  if isfield(cfg, 'parameter') && ~strcmp(cfg.parameter,'powspctrm')
    % freqdesctiptives will only work on the powspctrm field
    % hence a temporary copy of the data is needed
    tempdata.dimord    = data.dimord;
    tempdata.freq      = data.freq;
    tempdata.label     = data.label;
    tempdata.time      = data.time;
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

% Handle the bivariate case

% Check for bivariate metric with 'chan_chan' in the dimord
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;

% Check for bivariate metric with a labelcmb
haslabelcmb = isfield(data, 'labelcmb');

% check whether the bivariate metric is the one requested to plot
%shouldPlotCmb = (haslabelcmb && ...
%  size(data.(cfg.parameter),selchan(1)) == size(data.labelcmb,1)) ...
%  || isfull; % this should work because if dimord has multiple chans (so isfull=1)
%             % then we can never plot anything without reference channel
%             % this is different when haslabelcmb=1; then the parameter
%             % requested to plot might well be a simple powspctrm
%if (isfull || haslabelcmb) && shouldPlotCmb

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
      %data.(cfg.parameter) = reshape(mean(data.(cfg.parameter)(:,sel,:),2),[siz(1) 1 siz(3:end)]);
      sel1 = 1:siz(1);
      sel2 = sel;
      meandir = 2;
    elseif strcmp(cfg.directionality, 'outflow')
      %data.(cfg.parameter) = reshape(mean(data.(cfg.parameter)(sel,:,:),1),[siz(1) 1 siz(3:end)]);
      sel1 = sel;
      sel2 = 1:siz(1);
      meandir = 1;

    elseif strcmp(cfg.directionality, 'ff-fd')
      error('cfg.directionality = ''ff-fd'' is not supported anymore, you have to manually subtract the two before the call to ft_singleplotTFR');
    elseif strcmp(cfg.directionality, 'fd-ff')
      error('cfg.directionality = ''fd-ff'' is not supported anymore, you have to manually subtract the two before the call to ft_singleplotTFR');
    end %if directionality
  end %if ~isfull
end %handle the bivariate data

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

% Get physical x-axis range:
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

% Get physical y-axis range:
if strcmp(cfg.ylim,'maxmin')
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

% % test if X and Y are linearly spaced (to within 10^-12): % FROM UIMAGE
% x = data.(xparam)(xmin:xmax);
% y = data.(yparam)(ymin:ymax);
% dx = min(diff(x));  % smallest interval for X
% dy = min(diff(y));  % smallest interval for Y
% evenx = all(abs(diff(x)/dx-1)<1e-12);     % true if X is linearly spaced
% eveny = all(abs(diff(y)/dy-1)<1e-12);     % true if Y is linearly spaced
%
% % masking only possible for evenly spaced axis
% if strcmp(cfg.masknans, 'yes') && (~evenx || ~eveny)
%   warning('(one of the) axis are not evenly spaced -> nans cannot be masked out -> cfg.masknans is set to ''no'';')
%   cfg.masknans = 'no';
% end
%
% if ~isempty(cfg.maskparameter) && (~evenx || ~eveny)
%   warning('(one of the) axis are not evenly spaced -> no masking possible -> cfg.maskparameter cleared')
%   cfg.maskparameter = [];
% end

% perform channel selection
selchannel = ft_channelselection(cfg.channel, data.label);
sellab     = match_str(data.label, selchannel);

% cfg.maskparameter only possible for single channel
if length(sellab) > 1 && ~isempty(cfg.maskparameter)
  warning('no masking possible for average over multiple channels -> cfg.maskparameter cleared')
  cfg.maskparameter = [];
end

% get dimord dimensions
ydim = find(strcmp(yparam, dimtok));
xdim = find(strcmp(xparam, dimtok));
zdim = setdiff(1:length(dimtok), [ydim xdim]); % all other dimensions

% and permute
dat = data.(cfg.parameter);
dat = permute(dat, [zdim(:)' ydim xdim]);
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
    mask = squeeze(mask(sellab, ymin:ymax, xmin:xmax));
  elseif cfg.maskalpha == 1
    mask = squeeze(mask(sellab, ymin:ymax, xmin:xmax));
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
datamatrix = reshape(mean(dat, 1), [siz(2:end) 1]);
xvector    = data.(xparam)(xmin:xmax);
yvector    = data.(yparam)(ymin:ymax);

% Get physical z-axis range (color axis):
if strcmp(cfg.zlim,'maxmin')
  zmin = min(datamatrix(:));
  zmax = max(datamatrix(:));
elseif strcmp(cfg.zlim,'maxabs')
  zmin = -max(abs(datamatrix(:)));
  zmax = max(abs(datamatrix(:)));
elseif strcmp(cfg.zlim,'zeromax')
  zmin = 0;
  zmax = max(datamatrix(:));
elseif strcmp(cfg.zlim,'minzero')
  zmin = min(datamatrix(:));
  zmax = 0;
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% set colormap
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('singleplotTFR(): Colormap must be a n x 3 matrix'); end
  set(gcf,'colormap',cfg.colormap);
end

% Draw plot (and mask NaN's if requested):
cla
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
% set(gca,'Color','k')

if isequal(cfg.colorbar,'yes')
  % tag the colorbar so we know which axes are colorbars
  colorbar('tag', 'ft-colorbar');
end

% Set adjust color axis
if strcmp('yes',cfg.hotkeys)
  %  Attach data and cfg to figure and attach a key listener to the figure
  set(gcf, 'KeyPressFcn', {@key_sub, zmin, zmax})
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  % first, attach data to the figure with the current axis handle as a name
  dataname = fixname(num2str(double(gca)));
  setappdata(gcf,dataname,data);
  set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg}, 'event', 'WindowButtonMotionFcn'});
  %   set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
  %   set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
  %   set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'callback', {@select_topoplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
end

% Create title text containing channel name(s) and channel number(s):
if length(sellab) == 1
  t = [char(cfg.channel) ' / ' num2str(sellab) ];
else
  t = sprintf('mean(%0s)', join_str(',', cfg.channel));
end
h = title(t,'fontsize', cfg.fontsize);

% set the figure window title, add channel labels if number is small
if isempty(get(gcf, 'Name'))
  if length(sellab) < 5
    chans = join_str(',', cfg.channel);
  else
    chans = '<multiple channels>';
  end
  if isfield(cfg,'dataname')
    if iscell(cfg.dataname)
      dataname = cfg.dataname{1};
    else
      dataname = cfg.dataname;
    end
  elseif nargin > 1
    dataname = inputname(2);
  else % data provided through cfg.inputfile
    dataname = cfg.inputfile;
  end
  if isempty(cfg.figurename)
    set(gcf, 'Name', sprintf('%d: %s: %s (%s)', double(gcf), mfilename, dataname, chans));
    set(gcf, 'NumberTitle', 'off');
  else
    set(gcf, 'name', cfg.figurename);
    set(gcf, 'NumberTitle', 'off');
  end
end
axis tight;
hold off;

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance

% add a menu to the figure, but only if the current figure does not have subplots
% also, delete any possibly existing previous menu, this is safe because delete([]) does nothing
delete(findobj(gcf, 'type', 'uimenu', 'label', 'FieldTrip'));
if numel(findobj(gcf, 'type', 'axes', '-not', 'tag', 'ft-colorbar')) <= 1
  ftmenu = uimenu(gcf, 'Label', 'FieldTrip');
  uimenu(ftmenu, 'Label', 'Show pipeline',  'Callback', {@menu_pipeline, cfg});
  uimenu(ftmenu, 'Label', 'About',  'Callback', @menu_about);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotTFR(cfg, varargin)
% first to last callback-input of ft_select_range is range
% last callback-input of ft_select_range is contextmenu label, if used
range = varargin{end-1};
varargin = varargin(1:end-2); % remove range and last

% get appdata belonging to current axis
dataname = fixname(num2str(double(gca)));
data = getappdata(gcf, dataname);

if isfield(cfg, 'inputfile')
  % the reading has already been done and varargin contains the data
  cfg = rmfield(cfg, 'inputfile');
end

% make sure the topo displays all channels, not just the ones in this
% singleplot
cfg.channel = 'all';

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
ft_topoplotTFR(cfg, data);

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
