function [cfg] = ft_multiplotTFR(cfg, data)

% ft_multiplotTFR plots time-frequency representations of power or coherence in a 
% topographical layout. The plots of the indivual sensors are arranged according 
% to their location specified in the layout.
%
% Use as:
%   ft_multiplotTFR(cfg, data)
%
% The data can be a time-frequency representation of power or coherence that 
% was computed using the FT_FREQANALYSIS or FT_FREQDESCRIPTIVES functions.
%
% The configuration can have the following parameters:
% cfg.xparam           = field to be plotted on x-axis (default depends on data.dimord)
%                        'time'
% cfg.yparam           = field to be plotted on y-axis (default depends on data.dimord)
%                        'freq'
% cfg.zparam           = field to be represented as color (default depends on data.dimord)
%                        'powspctrm' or 'cohspctrm' 
% cfg.maskparameter    = field in the data to be used for opacity masking of data
% cfg.maskstyle        = style used to mask nans, 'opacity' or 'saturation' (default = 'opacity')
%                        use 'saturation' when saving to vector-format (like *.eps) to avoid all 
%                        sorts of image-problems (currently only possible with a white backgroud)
% cfg.xlim             = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim             = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.zlim             = 'maxmin','maxabs' or [zmin zmax] (default = 'maxmin')
% cfg.channel          = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
% cfg.cohrefchannel    = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline         = 'yes','no' or [time1 time2] (default = 'no'), see FT_FREQBASELINE
% cfg.baselinetype     = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials           = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.box              = 'yes', 'no' (default = 'no' if maskparameter given default = 'yes')
%                        Draw a box around each graph
% cfg.colorbar         = 'yes', 'no' (default = 'no')
% cfg.colormap         = any sized colormap, see COLORMAP
% cfg.comment          = string of text (default = date + zlimits)
%                        Add 'comment' to graph (according to COMNT in the layout)
% cfg.showlabels       = 'yes', 'no' (default = 'no')
% cfg.showoutline      = 'yes', 'no' (default = 'no')
% cfg.fontsize         = font size of comment and labels (if present) (default = 8)
% cfg.interactive      = Interactive plot 'yes' or 'no' (default = 'no')
%                        In a interactive plot you can select areas and produce a new
%                        interactive plot when a selected area is clicked. Multiple areas 
%                        can be selected by holding down the SHIFT key.
% cfg.masknans         = 'yes' or 'no' (default = 'yes')
% cfg.renderer         = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
% cfg.layout           = specify the channel layout for plotting using one of 
%                        the following ways:
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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
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
% cfg.xparam
% cfg.zparam

%
% This function depends on FT_FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype, documented

% Copyright (C) 2003-2006, Ole Jensen
% Copyright (C) 2007-2011, Roemer van der Meij & Jan-Mathijs Schoffelen
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

ft_defaults

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'unused', {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'zlim',  'absmax',  'maxabs'});

clf

% set default for inputfile
if ~isfield(cfg, 'inputfile'),      cfg.inputfile = [];                end

hasdata      = (nargin>1);
hasinputfile = ~isempty(cfg.inputfile);

if hasdata && hasinputfile
  error('cfg.inputfile should not be used in conjunction with giving input data to this function');
end

if hasdata
  % do nothing
elseif hasinputfile
  data = loadvar(cfg.inputfile, 'data');
end

% Set the defaults:
if ~isfield(cfg,'baseline'),        cfg.baseline = 'no';               end
if ~isfield(cfg,'baselinetype'),    cfg.baselinetype = 'absolute';     end
if ~isfield(cfg,'trials'),          cfg.trials = 'all';                end
if ~isfield(cfg,'xlim'),            cfg.xlim = 'maxmin';               end
if ~isfield(cfg,'ylim'),            cfg.ylim = 'maxmin';               end
if ~isfield(cfg,'zlim'),            cfg.zlim = 'maxmin';               end
if ~isfield(cfg,'colorbar'),        cfg.colorbar = 'no';               end
if ~isfield(cfg,'comment'),         cfg.comment = date;                end
if ~isfield(cfg,'showlabels'),      cfg.showlabels = 'no';             end
if ~isfield(cfg,'showoutline'),     cfg.showoutline = 'no';            end
if ~isfield(cfg,'channel'),         cfg.channel = 'all';               end
if ~isfield(cfg,'fontsize'),        cfg.fontsize = 8;                  end
if ~isfield(cfg,'interactive'),     cfg.interactive = 'no';            end
if ~isfield(cfg,'renderer'),        cfg.renderer = [];                 end % let matlab decide on default
if ~isfield(cfg,'masknans'),        cfg.masknans = 'yes';              end
if ~isfield(cfg,'maskparameter'),   cfg.maskparameter = [];            end
if ~isfield(cfg,'maskstyle'),       cfg.maskstyle = 'opacity';         end
if ~isfield(cfg,'matrixside'),      cfg.matrixside = 'feedforward';    end
if ~isfield(cfg,'channel'),         cfg.channel       = 'all';         end
if ~isfield(cfg,'box')             
  if ~isempty(cfg.maskparameter)
    cfg.box = 'yes';
  else
    cfg.box = 'no';
  end
end

% for backward compatibility with old data structures
data   = ft_checkdata(data, 'datatype', 'freq');
dimord = data.dimord;
dimtok = tokenize(dimord, '_');

% Set x/y/zparam defaults
if ~sum(ismember(dimtok, 'time'))
  error('input data needs a time dimension');
else
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
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
hasrpt = sum(ismember(dimtok, {'rpt' 'subj'}));
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
  if isfield(cfg, 'zparam') && ~strcmp(cfg.zparam,'powspctrm')
    % freqdesctiptives will only work on the powspctrm field
    % hence a temporary copy of the data is needed
    tempdata.dimord    = data.dimord;
    tempdata.freq      = data.freq;
    tempdata.label     = data.label;
    tempdata.powspctrm = data.(cfg.zparam);
    tempdata.cfg       = data.cfg;
    tempdata           = ft_freqdescriptives(tmpcfg, tempdata);
    data.(cfg.zparam)  = tempdata.powspctrm;
    clear tempdata
  else
    data = ft_freqdescriptives(tmpcfg, data);
  end
  dimord = data.dimord;
  dimtok = tokenize(dimord, '_');
end % if hasrpt

% Read or create the layout that will be used for plotting:
lay = ft_prepare_layout(cfg, data);
cfg.layout = lay;

% Apply baseline correction:
if ~strcmp(cfg.baseline, 'no')
  data = ft_freqbaseline(cfg, data);
end

% Handle the bivariate case

% Check for bivariate metric with 'chan_chan' in the dimord
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;

% Check for bivariate metric with a labelcmb
haslabelcmb = isfield(data, 'labelcmb');

if (isfull || haslabelcmb) && isfield(data, cfg.zparam)
  % A reference channel is required:
  if ~isfield(cfg, 'cohrefchannel')
    error('no reference channel is specified');
  end
  
  % check for cohrefchannel being part of selection
  if ~strcmp(cfg.cohrefchannel,'gui')
    if (isfull      && ~any(ismember(data.label, cfg.cohrefchannel))) || ...
       (haslabelcmb && ~any(ismember(data.labelcmb(:), cfg.cohrefchannel)))
      error('cfg.cohrefchannel is a not present in the (selected) channels)')
    end
  end
  
  % Interactively select the reference channel
  if strcmp(cfg.cohrefchannel, 'gui')
    % Open a single figure with the channel layout, the user can click on a reference channel
    h = clf;
    ft_plot_lay(lay, 'box', false);
    title('Select the reference channel by dragging a selection window, more than 1 channel can be selected...');
    % add the channel information to the figure
    info       = guidata(gcf);
    info.x     = lay.pos(:,1);
    info.y     = lay.pos(:,2);
    info.label = lay.label;
    guidata(h, info);
    %set(gcf, 'WindowButtonUpFcn', {@ft_select_channel, 'callback', {@select_topoplotER, cfg, data}});
    set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
    return
  end
  
  if ~isfull,
    % Convert 2-dimensional channel matrix to a single dimension:
    if isempty(cfg.matrixside)
      sel1 = strmatch(cfg.cohrefchannel, data.labelcmb(:,2));
      sel2 = strmatch(cfg.cohrefchannel, data.labelcmb(:,1));
    elseif strcmp(cfg.matrixside, 'feedforward')
      sel1 = [];
      sel2 = strmatch(cfg.cohrefchannel, data.labelcmb(:,1));
    elseif strcmp(cfg.matrixside, 'feedback')
      sel1 = strmatch(cfg.cohrefchannel, data.labelcmb(:,2));
      sel2 = [];
    end
    fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.zparam);
    data.(cfg.zparam) = data.(cfg.zparam)([sel1;sel2],:,:);
    data.label     = [data.labelcmb(sel1,1);data.labelcmb(sel2,2)];
    data.labelcmb  = data.labelcmb([sel1;sel2],:);
    data           = rmfield(data, 'labelcmb');
  else
    % General case
    sel               = match_str(data.label, cfg.cohrefchannel);
    siz               = [size(data.(cfg.zparam)) 1];
    if strcmp(cfg.matrixside, 'feedback') || isempty(cfg.matrixside)
      %FIXME the interpretation of 'feedback' and 'feedforward' depend on
      %the definition in the bivariate representation of the data
      %data.(cfg.zparam) = reshape(mean(data.(cfg.zparam)(:,sel,:),2),[siz(1) 1 siz(3:end)]);
      sel1 = 1:siz(1);
      sel2 = sel;
      meandir = 2;
    elseif strcmp(cfg.matrixside, 'feedforward')
      %data.(cfg.zparam) = reshape(mean(data.(cfg.zparam)(sel,:,:),1),[siz(1) 1 siz(3:end)]);
      sel1 = sel;
      sel2 = 1:siz(1);
      meandir = 1;

    elseif strcmp(cfg.matrixside, 'ff-fd')
      error('cfg.matrixside = ''ff-fd'' is not supported anymore, you have to manually subtract the two before the call to ft_topoplotER');
    elseif strcmp(cfg.matrixside, 'fd-ff')
      error('cfg.matrixside = ''fd-ff'' is not supported anymore, you have to manually subtract the two before the call to ft_topoplotER');
    end %if matrixside
  end %if ~isfull
end %handle the bivariate data


% Get physical x-axis range:
if strcmp(cfg.xlim,'maxmin')
  xmin = min(data.(cfg.xparam));
  xmax = max(data.(cfg.xparam));
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Replace value with the index of the nearest bin
if ~isempty(cfg.xparam)
  xmin = nearest(data.(cfg.xparam), xmin);
  xmax = nearest(data.(cfg.xparam), xmax);
end

% Get physical y-axis range:
if strcmp(cfg.ylim,'maxmin')
  ymin = min(data.(cfg.yparam));
  ymax = max(data.(cfg.yparam));
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Replace value with the index of the nearest bin
if ~isempty(cfg.yparam)
  ymin = nearest(data.(cfg.yparam), ymin);
  ymax = nearest(data.(cfg.yparam), ymax);
end

% test if X and Y are linearly spaced (to within 10^-12): % FROM UIMAGE
x = data.(cfg.xparam)(xmin:xmax);
y = data.(cfg.yparam)(ymin:ymax);
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

dat = data.(cfg.zparam);
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
  if isfull
    mask = mask(sel1, sel2, ymin:ymax, xmin:xmax);
    mask = nanmean(nanmean(nanmean(mask, meandir), 4), 3);
  elseif haslabelcmb
    mask = mask(sellab, ymin:ymax, xmin:xmax);
    mask = nanmean(nanmean(mask, 3), 2);
  else
    mask = mask(sellab, ymin:ymax, xmin:xmax);
    mask = nanmean(nanmean(mask, 3), 2);
  end
end

% Select the channels in the data that match with the layout:
[seldat, sellay] = match_str(label, lay.label);
if isempty(seldat)
  error('labels in data and labels in layout do not match'); 
end

datavector = dat(seldat,:,:);
if ~isempty(cfg.maskparameter)
  maskvector = mask(seldat,:,:);
end

% Select x and y coordinates and labels of the channels in the data
chanX = lay.pos(sellay, 1);
chanY = lay.pos(sellay, 2);
chanWidth  = lay.width(sellay);
chanHeight = lay.height(sellay);
chanLabels = lay.label(sellay);


% Get physical z-axis range (color axis):
if strcmp(cfg.zlim,'maxmin')
  zmin = min(datavector(:));
  zmax = max(datavector(:));
elseif strcmp(cfg.zlim,'maxabs')
  zmin = -max(abs(datavector(:)));
  zmax = max(abs(datavector(:)));
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

hold on;

if isfield(lay, 'outline') && strcmp(cfg.showoutline, 'yes')
  for i=1:length(lay.outline)
    if ~isempty(lay.outline{i})
      tmpX = lay.outline{i}(:,1);
      tmpY = lay.outline{i}(:,2);
      h = line(tmpX, tmpY);
      set(h, 'color', 'k');
      set(h, 'linewidth', 2);
    end
  end
end

% set colormap
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('multiplotTFR(): Colormap must be a n x 3 matrix'); end
  set(gcf,'colormap',cfg.colormap);
end;

% Plot channels:
for k=1:length(seldat)
  % Get cdata:
  cdata = squeeze(datavector(k,:,:));
  if ~isempty(cfg.maskparameter)
    mdata = squeeze(maskvector(k,:,:));
  end
  
  % Get axes for this panel
  xas = (chanX(k) + linspace(0,1,size(cdata,2))*chanWidth(k)) - chanWidth(k)/2;
  yas = (chanY(k) + linspace(0,1,size(cdata,1))*chanHeight(k)) - chanHeight(k)/2;
  
  % Draw plot (and mask Nan's with maskfield if requested)
  if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = double(mask);
    ft_plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = mask .* mdata;
    mask = double(mask);
    ft_plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  elseif isequal(cfg.masknans,'no') && ~isempty(cfg.maskparameter)
    mask = mdata;
    mask = double(mask);
    ft_plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  else
    ft_plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip')
  end
  % Currently the handle isn't being used below, this is here for possible use in the future
  h = findobj('tag','cip');
  
  % Draw box around plot
  if strcmp(cfg.box,'yes')
    xstep = xas(2) - xas(1); ystep = yas(2) - yas(1);
    xvalmin(1:length(yas)+2) = min(xas)-(0.5*xstep); xvalmax(1:length(yas)+2) = max(xas)+(0.5*xstep); yvalmin(1:length(xas)+2) = min(yas)-(0.5*ystep); yvalmax(1:length(xas)+2) = max(yas)+(0.5*ystep);
    xas2 = [xvalmin(1) xas xvalmax(1)]; yas2 = [yvalmin(1) yas yvalmax(1)];
    ft_plot_box([min(xas2) max(xas2) min(yas2) max(yas2)])
  end

  % Draw channel labels:
  if strcmp(cfg.showlabels,'yes')
    ft_plot_text(chanX(k)-chanWidth(k)/2, chanY(k)+chanHeight(k)/2, sprintf(' %0s\n ', chanLabels{k}), 'Fontsize', cfg.fontsize);
  end
end

% write comment:
k = cellstrmatch('COMNT',lay.label);
if ~isempty(k)
  comment = cfg.comment;
  comment = sprintf('%0s\nxlim=[%.3g %.3g]', comment, xmin, xmax);
  comment = sprintf('%0s\nylim=[%.3g %.3g]', comment, ymin, ymax);
  comment = sprintf('%0s\nzlim=[%.3g %.3g]', comment, zmin, zmax);
  ft_plot_text(lay.pos(k,1), lay.pos(k,2), sprintf(comment), 'Fontsize', cfg.fontsize);
end

% plot scale:
k = cellstrmatch('SCALE',lay.label);
if ~isempty(k)
  % Get average cdata across channels:
  cdata = squeeze(mean(datavector, 1));
  
  % Get axes for this panel:
  xas = (lay.pos(k,1) + linspace(0,1,size(cdata,2))*lay.width(k));
  yas = (lay.pos(k,2) + linspace(0,1,size(cdata,1))*lay.height(k));
  
  % Draw plot (and mask Nan's with maskfield if requested)
  if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = double(mask);
    ft_plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter)
    mask = ~isnan(cdata);
    mask = mask .* mdata;
    mask = double(mask);
    ft_plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  elseif isequal(cfg.masknans,'no') && ~isempty(cfg.maskparameter)
    mask = mdata;
    mask = double(mask);
    ft_plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
  else
    ft_plot_matrix(xas, yas, cdata,'clim',[zmin zmax],'tag','cip')
  end
  % Currently the handle isn't being used below, this is here for possible use in the future
  h = findobj('tag','cip');

end


% plot colorbar:
if isfield(cfg, 'colorbar') && (strcmp(cfg.colorbar, 'yes'))
  colorbar;
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
    % add the channel information to the figure
    info       = guidata(gcf);
    info.x     = lay.pos(:,1);
    info.y     = lay.pos(:,2);
    info.label = lay.label;
    guidata(gcf, info);

    set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, data}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, data}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotTFR, cfg, data}, 'event', 'WindowButtonMotionFcn'});
end

axis tight
axis off
if strcmp(cfg.box, 'yes')
  abc = axis;
  axis(abc + [-1 +1 -1 +1]*mean(abs(abc))/10)
end
orient landscape
hold off

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l = cellstrmatch(str,strlist)
l = [];
for k=1:length(strlist)
  if strcmp(char(str),char(strlist(k)))
    l = [l k];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called by ft_select_channel in case cfg.cohrefchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_multiplotTFR(label, cfg, varargin)
cfg.cohrefchannel = label;
fprintf('selected cfg.cohrefchannel = ''%s''\n', join(',', cfg.cohrefchannel));
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_multiplotTFR(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotTFR(label, cfg, varargin)
if ~isempty(label)
  cfg.xlim = 'maxmin';
  cfg.ylim = 'maxmin';
  cfg.channel = label;
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
