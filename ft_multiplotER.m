function [cfg] = ft_multiplotER(cfg, varargin)

% ft_multiplotER plots the event-related fields or potentials versus time
% or of oscillatory activity (power or coherence) versus frequency. Multiple
% datasets can be overlayed.  The plots are arranged according to their
% location specified in the layout.
%
% Use as:
%   ft_multiplotER(cfg, data)
%   ft_multiplotER(cfg, data, data2, ..., dataN)
%
% The data can be an ERP/ERF produced by FT_TIMELOCKANALYSIS, a powerspectrum
% produced by FT_FREQANALYSIS or a coherencespectrum produced by FT_FREQDESCRIPTIVES.
% If you specify multiple datasets they must contain the same channels, etc.
%
% The configuration can have the following parameters:
% cfg.xparam        = field to be plotted on x-axis (default depends on data.dimord)
%                     'time' or 'freq'
% cfg.zparam        = field to be plotted on y-axis (default depends on data.dimord)
%                     'avg', 'powspctrm' or 'cohspctrm'
% cfg.maskparameter = field in the first dataset to be used for marking significant data
% cfg.maskstyle     = style used for masking of data, 'box', 'thickness' or 'saturation' (default = 'box')
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
% cfg.cohrefchannel = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see FT_TIMELOCKBASELINE or FT_FREQBASELINE
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.axes          = 'yes', 'no' (default = 'yes')
%                     Draw x- and y-axes for each graph
% cfg.box           = 'yes', 'no' (default = 'no')
%                     Draw a box around each graph
% cfg.comment       = string of text (default = date + colors)
%                     Add 'comment' to graph (according to COMNT in the layout)
% cfg.showlabels    = 'yes', 'no' (default = 'no')
% cfg.showoutline   = 'yes', 'no' (default = 'no')
% cfg.fontsize      = font size of comment and labels (if present) (default = 8)
% cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'no')
%                     In a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. Multiple areas
%                     can be selected by holding down the SHIFT key.
% cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
% cfg.linestyle     = linestyle/marker type, see options of the matlab PLOT function (default = '-')
%                     can be a single style for all datasets, or a cell-array containing one style for each dataset
% cfg.linewidth     = linewidth in points (default = 0.5)
% cfg.graphcolor    = color(s) used for plotting the dataset(s) (default = 'brgkywrgbkywrgbkywrgbkyw')
%                     alternatively, colors can be specified as Nx3 matrix of RGB values
%
% cfg.layout        = specify the channel layout for plotting using one of
%                     the following ways:
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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure. For this particular function, the
% data should be provided as a cell array.
%
% See also:
%   FT_MULTIPLOTTFR, FT_SINGLEPLOTER, FT_SINGLEPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR,
%   FT_PREPARE_LAYOUT

% Undocumented local options:
% cfg.layoutname

%
% This function depends on FT_TIMELOCKBASELINE which has the following options:
% cfg.baseline, documented
% cfg.channel
% cfg.baselinewindow
% cfg.previous
% cfg.version
%
% This function depends on FT_FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype

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
cfg = ft_checkconfig(cfg, 'unused',  {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim', 'absmax', 'maxabs'});

cla

% set default for inputfile
if ~isfield(cfg, 'inputfile'),  cfg.inputfile = [];    end

hasdata      = nargin>1;
hasinputfile = ~isempty(cfg.inputfile);

if hasdata && hasinputfile
  error('cfg.inputfile should not be used in conjunction with giving input data to this function');
end

if hasdata
  % do nothing
elseif hasinputfile
  if ~ischar(cfg.inputfile)
    cfg.inputfile = {cfg.inputfile};
  end
  for i = 1:numel(cfg.inputfile)
    varargin{i} = loadvar(cfg.inputfile{i}, 'data'); % read datasets
  end
  if isfield(cfg, 'interactive') && strcmp(cfg.interactive, 'yes'),
    warning('switching off interactive mode, this is not supported when loading an inputfile from disk');
  end
end

% set the defaults:
if ~isfield(cfg,'baseline'),      cfg.baseline      = 'no';                        end
if ~isfield(cfg,'trials'),        cfg.trials        = 'all';                       end
if ~isfield(cfg,'xlim'),          cfg.xlim          = 'maxmin';                    end
if ~isfield(cfg,'ylim'),          cfg.ylim          = 'maxmin';                    end
if ~isfield(cfg,'comment'),       cfg.comment       = strcat([date '\n']);         end
if ~isfield(cfg,'axes'),          cfg.axes          = 'yes';                       end
if ~isfield(cfg,'showlabels'),    cfg.showlabels    = 'no';                        end
if ~isfield(cfg,'showoutline'),   cfg.showoutline   = 'no';                        end
if ~isfield(cfg,'box'),           cfg.box           = 'no';                        end
if ~isfield(cfg,'fontsize'),      cfg.fontsize      = 8;                           end
if ~isfield(cfg,'graphcolor'),    cfg.graphcolor    = 'brgkywrgbkywrgbkywrgbkyw';  end
if ~isfield(cfg,'interactive'),   cfg.interactive   = 'no';                        end
if ~isfield(cfg,'renderer'),      cfg.renderer      = [];                          end
if ~isfield(cfg,'maskparameter'), cfg.maskparameter = [];                          end
if ~isfield(cfg,'linestyle'),     cfg.linestyle     = '-';                         end
if ~isfield(cfg,'linewidth'),     cfg.linewidth     = 0.5;                         end
if ~isfield(cfg,'maskstyle'),     cfg.maskstyle     = 'box';                       end
if ~isfield(cfg,'channel'),       cfg.channel       = 'all';                       end
if ~isfield(cfg, 'matrixside'),   cfg.matrixside    = '';                          end

Ndata = numel(varargin);

%FIXME rename matrixside and cohrefchannel in more meaningful options
if ischar(cfg.graphcolor)
  GRAPHCOLOR = ['k' cfg.graphcolor];
elseif isnumeric(cfg.graphcolor)
  GRAPHCOLOR = [0 0 0; cfg.graphcolor];
end

% check for linestyle being a cell-array, check it's length, and lengthen it if does not have enough styles in it
if ischar(cfg.linestyle)
  cfg.linestyle = {cfg.linestyle};
end

if Ndata  > 1
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

% ensure that the input is correct, also backward compatibility with old data structures:
for i=1:Ndata
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'timelock', 'freq'});
  dtype{i}    = ft_datatype(varargin{i});
  
  % this is needed for correct treatment of GRAPHCOLOR later on
  if nargin>1,
    if ~isempty(inputname(i+1))
      iname{i+1} = inputname(i+1);
    else
      iname{i+1} = ['input',num2str(i,'%02d')];
    end
  else
    iname{i+1} = cfg.inputfile{i};
  end
end

if Ndata>1,
  if ~all(strcmp(dtype{1}, dtype))
    error('input data are of different type; this is not supported');
  end
end
dtype  = dtype{1};
dimord = varargin{1}.dimord;
dimtok = tokenize(dimord, '_');

% Set x/y/zparam defaults according to datatype and dimord
switch dtype
  case 'timelock'
    if ~isfield(cfg, 'xparam'),      cfg.xparam = 'time';         end
    if ~isfield(cfg, 'yparam'),      cfg.yparam = '';             end
    if ~isfield(cfg, 'zparam'),      cfg.zparam = 'avg';          end
  case 'freq'
    if sum(ismember(dimtok, 'time'))
      if ~isfield(cfg, 'xparam'),    cfg.xparam = 'time';         end
      if ~isfield(cfg, 'yparam'),    cfg.yparam = 'freq';         end %FIXME
      if ~isfield(cfg, 'zparam'),    cfg.zparam = 'powspctrm';    end
    else
      if ~isfield(cfg, 'xparam'),    cfg.xparam = 'freq';         end
      if ~isfield(cfg, 'yparam'),    cfg.yparam = '';             end
      if ~isfield(cfg, 'zparam'),    cfg.zparam = 'powspctrm';    end
    end
  case 'comp'
    % not supported
  otherwise
    % not supported
end

% user specified own fields, but no yparam (which is not asked in help)
if isfield(cfg, 'xparam') && isfield(cfg, 'zparam') && ~isfield(cfg, 'yparam')
  cfg.yparam = '';
end

if isfield(cfg, 'channel') && isfield(varargin{1}, 'label')
  cfg.channel = ft_channelselection(cfg.channel, varargin{1}.label);
elseif isfield(cfg, 'channel') && isfield(varargin{1}, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(varargin{1}.labelcmb(:)));
end

% perform channel selection but only allow this when cfg.interactive = 'no'
if isfield(varargin{1}, 'label') && strcmp(cfg.interactive, 'no')
  selchannel = ft_channelselection(cfg.channel, varargin{1}.label);
elseif isfield(varargin{1}, 'labelcmb') && strcmp(cfg.interactive, 'no')
  selchannel = ft_channelselection(cfg.channel, unique(varargin{1}.labelcmb(:)));
end

% check whether rpt/subj is present and remove if necessary and whether
hasrpt = sum(ismember(dimtok, {'rpt' 'subj'}));
if strcmp(dtype, 'timelock') && hasrpt,
  tmpcfg        = [];
  tmpcfg.trials = cfg.trials;
  for i=1:Ndata
    varargin{i} = ft_timelockanalysis(tmpcfg, varargin{i});
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
    if isfield(cfg, 'zparam') && ~strcmp(cfg.zparam,'powspctrm')
      % freqdesctiptives will only work on the powspctrm field
      % hence a temporary copy of the data is needed
      tempdata.dimord    = varargin{i}.dimord;
      tempdata.freq      = varargin{i}.freq;
      tempdata.label     = varargin{i}.label;
      tempdata.powspctrm = varargin{i}.(cfg.zparam);
      tempdata.cfg       = varargin{i}.cfg;
      tempdata           = ft_freqdescriptives(tmpcfg, tempdata);
      varargin{i}.(cfg.zparam)  = tempdata.powspctrm;
      clear tempdata
    else
      varargin{i} = ft_freqdescriptives(tmpcfg, varargin{i});
    end
  end
  dimord = varargin{1}.dimord;
  dimtok = tokenize(dimord, '_');
end

% Read or create the layout that will be used for plotting
lay = ft_prepare_layout(cfg, varargin{1});
cfg.layout = lay;

% Apply baseline correction
if ~strcmp(cfg.baseline, 'no')
  for i=1:Ndata
    if strcmp(dtype, 'timelock') && strcmp(cfg.xparam, 'time')
      varargin{i} = ft_timelockbaseline(cfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(cfg.xparam, 'time')
      varargin{i} = ft_freqbaseline(cfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(cfg.xparam, 'freq')
      error('Baseline correction is not supported for spectra without a time dimension');
    else
      warning('Baseline correction not applied, please set cfg.xparam');
    end
  end
end

% Handle the bivariate case

% Check for bivariate metric with 'chan_chan' in the dimord
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;

% Check for bivariate metric with a labelcmb
haslabelcmb = isfield(varargin{1}, 'labelcmb');

if (isfull || haslabelcmb) && isfield(varargin{1}, cfg.zparam)
  % A reference channel is required:
  if ~isfield(cfg, 'cohrefchannel')
    error('no reference channel is specified');
  end
  
  % check for cohrefchannel being part of selection
  if ~strcmp(cfg.cohrefchannel,'gui')
    if (isfull      && ~any(ismember(varargin{1}.label, cfg.cohrefchannel))) || ...
        (haslabelcmb && ~any(ismember(varargin{1}.labelcmb(:), cfg.cohrefchannel)))
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
    set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotER, cfg, varargin{1}}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotER, cfg, varargin{1}}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_multiplotER, cfg, varargin{1}}, 'event', 'WindowButtonMotionFcn'});
    return
  end
  
  for i=1:Ndata
    if ~isfull,
      % Convert 2-dimensional channel matrix to a single dimension:
      if isempty(cfg.matrixside)
        sel1 = strmatch(cfg.cohrefchannel, varargin{i}.labelcmb(:,2), 'exact');
        sel2 = strmatch(cfg.cohrefchannel, varargin{i}.labelcmb(:,1), 'exact');
      elseif strcmp(cfg.matrixside, 'feedforward')
        sel1 = [];
        sel2 = strmatch(cfg.cohrefchannel, varargin{i}.labelcmb(:,1), 'exact');
      elseif strcmp(cfg.matrixside, 'feedback')
        sel1 = strmatch(cfg.cohrefchannel, varargin{i}.labelcmb(:,2), 'exact');
        sel2 = [];
      end
      fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.zparam);
      varargin{i}.(cfg.zparam) = varargin{i}.(cfg.zparam)([sel1;sel2],:,:);
      varargin{i}.label     = [varargin{i}.labelcmb(sel1,1);varargin{i}.labelcmb(sel2,2)];
      varargin{i}.labelcmb  = varargin{i}.labelcmb([sel1;sel2],:);
      varargin{i}           = rmfield(varargin{i}, 'labelcmb');
    else
      % General case
      sel               = match_str(varargin{i}.label, cfg.cohrefchannel);
      siz               = [size(varargin{i}.(cfg.zparam)) 1];
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
  end %for i
end %handle the bivariate data

% Get physical min/max range of x
if strcmp(cfg.xlim,'maxmin')
  % Find maxmin throughout all varargins:
  xmin = [];
  xmax = [];
  for i=1:length(varargin)
    xmin = min([xmin varargin{i}.(cfg.xparam)]);
    xmax = max([xmax varargin{i}.(cfg.xparam)]);
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Get the index of the nearest bin
for i=1:Ndata
  xidmin(i,1) = nearest(varargin{i}.(cfg.xparam), xmin);
  xidmax(i,1) = nearest(varargin{i}.(cfg.xparam), xmax);
end

%FIXME do something with yparam here
%technically should not be defined for multiplotER, but can be defined (and
%use ft_selectdata to average across frequencies

% Get physical y-axis range (ylim / zparam):
if strcmp(cfg.ylim,'maxmin')
  % Find maxmin throughout all varargins:
  ymin = [];
  ymax = [];
  for i=1:length(varargin)
    % Select the channels in the data that match with the layout:
    dat = [];
    dat = varargin{i}.(cfg.zparam);
    [seldat, sellay] = match_str(varargin{i}.label, lay.label);
    if isempty(seldat)
      error('labels in data and labels in layout do not match');
    end
    data = dat(seldat,:);
    ymin = min([ymin min(min(min(data)))]);
    ymax = max([ymax max(max(max(data)))]);
  end
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% convert the layout to Ole's style of variable names
X      = lay.pos(:,1);
Y      = lay.pos(:,2);
width  = lay.width;
height = lay.height;
Lbl    = lay.label;

% Create empty channel coordinates and labels arrays:
chanX(1:length(Lbl)) = NaN;
chanY(1:length(Lbl)) = NaN;
chanLabels = cell(1,length(Lbl));

hold on;
colorLabels = [];

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

% Plot each data set:
for i=1:Ndata
  % Make vector dat with one value for each channel
  dat    = varargin{i}.(cfg.zparam);
  xparam = varargin{i}.(cfg.xparam);
  
  % Take subselection of channels, this only works
  % in the interactive mode
  if exist('selchannel', 'var')
    sellab = match_str(varargin{i}.label, selchannel);
    label  = varargin{i}.label(sellab);
  else
    sellab = 1:numel(varargin{i}.label);
    label  = varargin{i}.label;
  end
  
  if ~isempty(cfg.yparam)
    if isfull
      dat = dat(sel1, sel2, ymin:ymax, xidmin(i):xidmax(i));
      dat = nanmean(nanmean(dat, meandir), 3);
      siz = size(dat);
      %FIXMEdat = reshape(dat, [siz(1:2) siz(4)]);
    elseif haslabelcmb
      dat = dat(sellab, ymin:ymax, xidmin(i):xidmax(i));
      dat = nanmean(dat, 2);
      siz = size(dat);
      dat = reshape(dat, [siz(1) siz(3)]);
    else
      dat = dat(sellab, ymin:ymax, xidmin(i):xidmax(i));
      dat = nanmean(nanmean(dat, 3), 2);
      siz = size(dat);
      dat = reshape(dat, [siz(1) siz(3)]);
    end
  else
    if isfull
      dat = dat(sel1, sel2, xidmin(i):xidmax(i));
      dat = nanmean(dat, meandir);
      %FIXME
    elseif haslabelcmb
      dat = dat(sellab, xidmin(i):xidmax(i));
    else
      dat = dat(sellab, xidmin(i):xidmax(i));
    end
  end
  xparam = xparam(xidmin(i):xidmax(i));
  
  % Select the channels in the data that match with the layout:
  [seldat, sellay] = match_str(label, cfg.layout.label);
  if isempty(seldat)
    error('labels in data and labels in layout do not match');
  end
  
  datamatrix = dat(seldat, :);
  
  % Select x and y coordinates and labels of the channels in the data
  layX = cfg.layout.pos(sellay,1);
  layY = cfg.layout.pos(sellay,2);
  layLabels = cfg.layout.label(sellay);
  
  if ~isempty(cfg.maskparameter)
    % one value for each channel, or one value for each channel-time point
    maskmatrix = varargin{1}.(cfg.maskparameter)(seldat,:);
    maskmatrix = maskmatrix(:,xidmin:xidmax);
  else
    % create an Nx0 matrix
    maskmatrix = zeros(length(seldat), 0);
  end
  
  if Ndata > 1
    if ischar(GRAPHCOLOR);        colorLabels = [colorLabels iname{i+1} '=' GRAPHCOLOR(i+1) '\n'];
    elseif isnumeric(GRAPHCOLOR); colorLabels = [colorLabels iname{i+1} '=' num2str(GRAPHCOLOR(i+1,:)) '\n'];
    end
  end
  
  if ischar(GRAPHCOLOR);        color = GRAPHCOLOR(i+1);
  elseif isnumeric(GRAPHCOLOR); color = GRAPHCOLOR(i+1,:);
  end
  
  for m=1:length(layLabels)
    % Plot ER
    plotWnd(xparam, datamatrix(m,:),[xmin xmax],[ymin ymax], layX(m), layY(m), width(m), height(m), layLabels(m), cfg, color, cfg.linestyle{i}, maskmatrix(m,:),i); %FIXME shouldn't this be replaced with a call to ft_plot_vector?
    
    if i==1,
      % Keep ER plot coordinates (at centre of ER plot), and channel labels (will be stored in the figure's UserData struct):
      chanX(m) = X(m) + 0.5 * width(m);
      chanY(m) = Y(m) + 0.5 * height(m);
      chanLabels{m} = Lbl{m};
    end
  end %for m
end %for i

% Add the colors of the different datasets to the comment:
cfg.comment = [cfg.comment colorLabels];

% Write comment text:
l = cellstrmatch('COMNT',Lbl);
if ~isempty(l)
  ft_plot_text(X(l),Y(l),sprintf(cfg.comment),'Fontsize',cfg.fontsize);
end

% Plot scales:
l = cellstrmatch('SCALE',Lbl);
if ~isempty(l)
  plotScales([xmin xmax],[ymin ymax],X(l),Y(l),width(1),height(1),cfg)
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  
  % add the channel information to the figure
  info       = guidata(gcf);
  info.x     = lay.pos(:,1);
  info.y     = lay.pos(:,2);
  info.label = lay.label;
  guidata(gcf, info);
  
  set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER, cfg, varargin{:}}, 'event', 'WindowButtonMotionFcn'});
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
function plotScales(xlim,ylim,xpos,ypos,width,height,cfg)
x1 =  xpos;
x2 =  xpos+width;
y1 =  ypos;
y2 =  ypos+width;
ft_plot_box([xpos xpos+width ypos ypos+height],'edgecolor','b')

if xlim(1) <=  0 && xlim(2) >= 0
  xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','b');
end

if ylim(1) <= 0 && ylim(2) >= 0
  xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','b');
end

ft_plot_text( x1,y1,num2str(xlim(1),3),'rotation',90,'HorizontalAlignment','Right','VerticalAlignment','middle','Fontsize',cfg.fontsize);
ft_plot_text( x2,y1,num2str(xlim(2),3),'rotation',90,'HorizontalAlignment','Right','VerticalAlignment','middle','Fontsize',cfg.fontsize);
ft_plot_text( x2,y1,num2str(ylim(1),3),'HorizontalAlignment','Left','VerticalAlignment','bottom','Fontsize',cfg.fontsize);
ft_plot_text( x2,y2,num2str(ylim(2),3),'HorizontalAlignment','Left','VerticalAlignment','bottom','Fontsize',cfg.fontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotWnd(x,y,xlim,ylim,xpos,ypos,width,height,label,cfg,color,style,mask,i)

% Clip out of bounds y values:
y(y > ylim(2)) = ylim(2);
y(y < ylim(1)) = ylim(1);

xs = xpos+width*(x-xlim(1))/(xlim(2)-xlim(1));
ys = ypos+height*(y-ylim(1))/(ylim(2)-ylim(1));

% Add boxes when masktyle is box, ft_plot_vector doesnt support boxes higher than ydata yet, so this code is left here
if i<2 && ~isempty(mask) && strcmp(cfg.maskstyle, 'box') % i stops box from being plotted more than once
  % determine how many boxes
  mask = mask(:)';
  mask = mask~=0;
  mask = diff([0 mask 0]);
  boxbeg = find(mask== 1);
  boxend = find(mask==-1)-1;
  
  numbox = length(boxbeg);
  for i = 1:numbox
    xmaskmin = xpos+width*(x(boxbeg(i))-xlim(1))/(xlim(2)-xlim(1));
    xmaskmax = xpos+width*(x(boxend(i))-xlim(1))/(xlim(2)-xlim(1));
    %plot([xmaskmin xmaskmax xmaskmax xmaskmin xmaskmin],[ypos ypos ypos+height ypos+height ypos],'r');
    hs = patch([xmaskmin xmaskmax xmaskmax xmaskmin xmaskmin],[ypos ypos ypos+height ypos+height ypos], [.6 .6 .6]);
    set(hs, 'EdgeColor', 'none');
  end
end

if isempty(mask) || (~isempty(mask) && strcmp(cfg.maskstyle,'box'))
  ft_plot_vector(xs, ys, 'color', color, 'style', style, 'linewidth', cfg.linewidth);
elseif ~isempty(mask) && ~strcmp(cfg.maskstyle,'box') % ft_plot_vector does not support boxes higher than ydata yet, so a separate option remains below
  ft_plot_vector(xs, ys, 'color', color, 'style', style, 'linewidth', cfg.linewidth, 'highlight', mask, 'highlightstyle', cfg.maskstyle);
end

if strcmp(cfg.showlabels,'yes')
  ft_plot_text(xpos,ypos+1.0*height,label,'Fontsize',cfg.fontsize);
end

% Draw axes:
if strcmp(cfg.axes,'yes') || strcmp(cfg.axes, 'xy')
  % Draw y axis
  xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','k');
  % Draw x axis
  xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','k');
  
elseif strcmp(cfg.axes,'x')
  % Draw x axis
  xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','k');
  
elseif strcmp(cfg.axes,'y')
  % Draw y axis
  xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
  ft_plot_vector(xs,ys,'color','k');
end

% Draw box around plot:
if strcmp(cfg.box,'yes')
  ft_plot_box([xpos xpos+width ypos ypos+height],'edgecolor','k');
end

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
% SUBFUNCTION which is called after selecting channels in case of cfg.cohrefchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_multiplotER(label, cfg, varargin)
if iscell(label)
  label = label{1};
end
cfg.cohrefchannel = label; %FIXME this only works with label being a string
fprintf('selected cfg.cohrefchannel = ''%s''\n', cfg.cohrefchannel);
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_multiplotER(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label, cfg, varargin)
if ~isempty(label)
  %cfg.xlim = 'maxmin';
  cfg.channel = label;
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
