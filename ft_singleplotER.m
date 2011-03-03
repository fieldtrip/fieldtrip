function [cfg] = ft_singleplotER(cfg, varargin)

% ft_singleplotER plots the event-related fields or potentials of a single channel
% or the average over multiple channels. Multiple datasets can be overlayed.
%
% Use as:
%   ft_sinlgeplotER(cfg, data)
%   ft_singleplotER(cfg, data1, data2, ..., dataN)
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
% cfg.maskparameter = field in the first dataset to be used for masking of data 
%                     (not possible for mean over multiple channels, or when input contains multiple subjects
%                     or trials)
% cfg.maskstyle     = style used for masking of data, 'box', 'thickness' or 'saturation' (default = 'box')
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
% cfg.cohrefchannel = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see FT_TIMELOCKBASELINE 
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.fontsize      = font size of title (default = 8)
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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also:
%   FT_SINGLEPLOTTFR, FT_MULTIPLOTER, FT_MULTIPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR

% This function depends on FT_TIMELOCKBASELINE which has the following options:
% cfg.baseline, documented
% cfg.channel
% cfg.baselinewindow
% cfg.previous
% cfg.version

% Copyright (C) 2003-2006, Ole Jensen
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
cfg = ft_checkconfig(cfg, 'renamed', {'channelindex',  'channel'});
cfg = ft_checkconfig(cfg, 'renamed', {'channelname',   'channel'});

cla

% set default for inputfile
cfg.inputfile = ft_getopt(cfg, 'inputfile', []);

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
cfg.baseline      = ft_getopt(cfg, 'baseline',    'no');
cfg.trials        = ft_getopt(cfg, 'trials',      'all');
cfg.xlim          = ft_getopt(cfg, 'xlim',        'maxmin'); 
cfg.ylim          = ft_getopt(cfg, 'ylim',        'maxmin');
cfg.comment       = ft_getopt(cfg, 'comment',     strcat([date '\n']));
cfg.axes          = ft_getopt(cfg,' axes',        'yes');
cfg.fontsize      = ft_getopt(cfg, 'fontsize',    8);
cfg.graphcolor    = ft_getopt(cfg, 'graphcolor',  'brgkywrgbkywrgbkywrgbkyw');
cfg.interactive   = ft_getopt(cfg, 'interactive',  'no');
cfg.renderer      = ft_getopt(cfg, 'renderer',     []);
cfg.maskparameter = ft_getopt(cfg, 'maskparameter',[]);
cfg.linestyle     = ft_getopt(cfg, 'linestyle',    '-');
cfg.linewidth     = ft_getopt(cfg, 'linewidth',    0.5);
cfg.maskstyle     = ft_getopt(cfg, 'maskstyle',    'box');
cfg.channel       = ft_getopt(cfg, 'channel',      'all');
cfg.matrixside    = ft_getopt(cfg, 'matrixside',   '');

Ndata = numel(varargin);

% interactive plotting is not allowed with more than 1 input
if Ndata >1 && strcmp(cfg.interactive, 'yes')
  error('interactive plotting is not supported with more than 1 input data set');
end

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

% ensure that the input is correct, also backward compatibility with old data structures:
dtype = cell(Ndata , 1);
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

if Ndata >1,
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
    cfg.xparam = ft_getopt(cfg,  'xparam', 'time'); 
    cfg.yparam = ft_getopt(cfg,  'yparam', '');     
    cfg.zparam = ft_getopt(cfg,  'zparam', 'avg');
  case 'freq'
    if sum(ismember(dimtok, 'time'))
      cfg.xparam = ft_getopt(cfg,  'xparam', 'time');
      cfg.yparam = ft_getopt(cfg,  'yparam', 'freq');%FIXME
      cfg.zparam = ft_getopt(cfg,  'zparam', 'powspctrm');
    else
      cfg.xparam = ft_getopt(cfg,  'xparam', 'freq'); 
      cfg.yparam = ft_getopt(cfg,  'yparam', '');     
      cfg.zparam = ft_getopt(cfg,  'zparam', 'powspctrm');
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
    error('cfg.cohrefchannel = ''gui'' is not supported in ft_singleplotER');
  end
  
  for i=1:Ndata 
    if ~isfull,
      % Convert 2-dimensional channel matrix to a single dimension:
      if isempty(cfg.matrixside)
        sel1 = strmatch(cfg.cohrefchannel, varargin{i}.labelcmb(:,2));
        sel2 = strmatch(cfg.cohrefchannel, varargin{i}.labelcmb(:,1));
      elseif strcmp(cfg.matrixside, 'feedforward')
        sel1 = [];
        sel2 = strmatch(cfg.cohrefchannel, varargin{i}.labelcmb(:,1));
      elseif strcmp(cfg.matrixside, 'feedback')
        sel1 = strmatch(cfg.cohrefchannel, varargin{i}.labelcmb(:,2));
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
  for i=1:Ndata 
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

hold on;
colorLabels = [];

% Plot each data set:
for i=1:Ndata 
  if isfield(varargin{1}, 'label')
    selchannel = ft_channelselection(cfg.channel, varargin{i}.label);
  elseif isfield(varargin{1}, 'labelcmb')
    selchannel = ft_channelselection(cfg.channel, unique(varargin{i}.labelcmb(:)));
  else
    error('the input data does not contain a label or labelcmb-field');
  end
  
  % Make vector dat with one value for each channel
  dat    = varargin{i}.(cfg.zparam);
  xparam = varargin{i}.(cfg.xparam);
  
  % Take subselection of channels
  % this works for bivariate data with labelcmb because at this point the
  % data has a label-field
  sellab = match_str(varargin{i}.label, selchannel);
  
  if ~isempty(cfg.yparam)
    if isfull
      dat = dat(sel1, sel2, ymin:ymax, xidmin(i):xidmax(i));
      dat = nanmean(nanmean(dat, meandir), 3);
      siz = size(dat);
      %FIXMEdat = reshape(dat, [siz(1:2) siz(4)]);
      dat = reshape(dat, [siz(1) siz(3)]);
      dat = dat(sellab, :);
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
      siz = size(dat);
      dat = reshape(dat, [siz(1) siz(3)]);
      dat = dat(sellab, :);
    elseif haslabelcmb
      dat = dat(sellab, xidmin(i):xidmax(i));
    else
      dat = dat(sellab, xidmin(i):xidmax(i));
    end
  end
  xvector    = xparam(xidmin(i):xidmax(i));
  datavector = reshape(mean(dat, 1), [1 numel(xvector)]); % average over channels

  % make mask
  if ~isempty(cfg.maskparameter)
    datmask = varargin{1}.(cfg.maskparameter);
    
    maskdatavector = reshape(mean(datmask(sellab,:),1), [1 numel(xvector)]);
  else
    maskdatavector = [];
  end

  if Ndata  > 1
    if ischar(GRAPHCOLOR);        colorLabels = [colorLabels iname{i+1} '=' GRAPHCOLOR(i+1) '\n'];
    elseif isnumeric(GRAPHCOLOR); colorLabels = [colorLabels iname{i+1} '=' num2str(GRAPHCOLOR(i+1,:)) '\n'];
    end
  end
  
  if ischar(GRAPHCOLOR);        color = GRAPHCOLOR(i+1);
  elseif isnumeric(GRAPHCOLOR); color = GRAPHCOLOR(i+1,:);
  end
    
  % Update ymin and ymax for the current data set:
  if strcmp(cfg.ylim,'maxmin')
    % Find maxmin for all varargins:
    ymin = zeros(1, Ndata );
    ymax = zeros(1, Ndata );
    for i=1:Ndata 
      % Select the channels in the data that match with the layout:
      ymin(i) = min(datavector);
      ymax(i) = max(datavector);
    end
    ymin = min(ymin);
    ymax = max(ymax);
  else
    ymin = cfg.ylim(1);
    ymax = cfg.ylim(2);
  end
 

  % only plot the mask once, for the first line (it's the same anyway for
  % all lines, and if plotted multiple times, it will overlay the others
  if i>1 && strcmp(cfg.maskstyle, 'box')
    ft_plot_vector(xvector, datavector, 'style', cfg.linestyle{i}, 'color', color, ...
   'linewidth', cfg.linewidth, 'hlim', cfg.xlim, 'vlim', cfg.ylim);  
  else
   ft_plot_vector(xvector, datavector, 'style', cfg.linestyle{i}, 'color', color, ...
   'highlight', maskdatavector, 'highlightstyle', cfg.maskstyle, 'linewidth', cfg.linewidth, ...
   'hlim', cfg.xlim, 'vlim', cfg.ylim);  
  end
end

% Set xlim and ylim:
xlim([xmin xmax]);
ylim([ymin ymax]);

% Make the figure interactive
if strcmp(cfg.interactive, 'yes')
  set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'WindowButtonMotionFcn'});
end

% Create title text containing channel name(s) and channel number(s):
if length(sellab) == 1
  t = [char(cfg.channel) ' / ' num2str(sellab) ];
else
  t = sprintf('mean(%0s)', join(',', cfg.channel));
end
h = title(t,'fontsize', cfg.fontsize);

% Set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called by ft_select_channel in case cfg.cohrefchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label, cfg, varargin)
cfg.cohrefchannel = label;
fprintf('selected cfg.cohrefchannel = ''%s''\n', cfg.cohrefchannel);
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_singleplotER(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotER(range, cfg, varargin)
cfg.comment = 'auto';
cfg.yparam = [];
cfg.xlim = range(1:2);
fprintf('selected cfg.xlim = [%f %f]\n', cfg.xlim(1), cfg.xlim(2));
p = get(gcf, 'Position');
f = figure;
set(f, 'Position', p);
ft_topoplotER(cfg, varargin{:});
