function [cfg] = ft_singleplotER(cfg, varargin)

% ft_singleplotER plots the event-related fields or potentials of a single channel
% or the average over multiple channels. multiple datasets can be overlayed.
%
% use as:
%   ft_singleplotER(cfg, data)
%   ft_singleplotER(cfg, data1, data2, ..., datan)
%
% the data can be an erp/erf produced by ft_timelockanalysis, a powerspectrum
% produced by ft_freqanalysis or a coherencespectrum produced by ft_freqdescriptives.
% if you specify multiple datasets they must contain the same channels, etc.
%
% the configuration can have the following parameters:
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
% cfg.channel       = nx1 cell-array with selection of channels (default = 'all'),
%                     see ft_channelselection for details
% cfg.refchannel    = name of reference channel for visualising connectivity, can be 'gui'
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see ft_timelockbaseline
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xn vector (default = 'all')
% cfg.fontsize      = font size of title (default = 8)
% cfg.hotkeys       = enables hotkeys (up/down/left/right arrows) for dynamic x/y axis translation (Ctrl+) and zoom adjustment
% cfg.interactive   = interactive plot 'yes' or 'no' (default = 'no')
%                     in a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. multiple areas
%                     can be selected by holding down the shift key.
% cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
% cfg.linestyle     = linestyle/marker type, see options of the matlab plot function (default = '-')
%                     can be a single style for all datasets, or a cell-array containing one style for each dataset
% cfg.linewidth     = linewidth in points (default = 0.5)
% cfg.graphcolor    = color(s) used for plotting the dataset(s) (default = 'brgkywrgbkywrgbkywrgbkyw')
%                     alternatively, colors can be specified as nx3 matrix of rgb values
%
% to facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% if you specify this option the input data will be read from a *.mat
% file on disk. this mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% see also:
%   ft_singleplotTFR, ft_multiplotER, ft_multiplotTFR, ft_topoplotER, ft_topoplotTFR

% this function depends on ft_timelockbaseline which has the following options:
% cfg.baseline, documented
% cfg.channel
% cfg.baselinewindow
% cfg.previous
% cfg.version

% copyright (c) 2003-2006, ole jensen
%
% this file is part of fieldtrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    fieldtrip is free software: you can redistribute it and/or modify
%    it under the terms of the gnu general public license as published by
%    the free software foundation, either version 3 of the license, or
%    (at your option) any later version.
%
%    fieldtrip is distributed in the hope that it will be useful,
%    but without any warranty; without even the implied warranty of
%    merchantability or fitness for a particular purpose.  see the
%    gnu general public license for more details.
%
%    you should have received a copy of the gnu general public license
%    along with fieldtrip. if not, see <http://www.gnu.org/licenses/>.
%
% $id: ft_singleplotER.m 3147 2011-03-17 12:38:09z jansch $

ft_defaults

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'unused',  {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim', 'absmax', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'matrixside',   'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'matrixside',   'feedback',    'inflow'});
cfg = ft_checkconfig(cfg, 'renamed', {'channelindex',  'channel'});
cfg = ft_checkconfig(cfg, 'renamed', {'channelname',   'channel'});
cfg = ft_checkconfig(cfg, 'renamed', {'cohrefchannel', 'refchannel'});


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
cfg.hotkeys       = ft_getopt(cfg, 'hotkeys', 'no');
cfg.interactive   = ft_getopt(cfg, 'interactive',  'no');
cfg.renderer      = ft_getopt(cfg, 'renderer',     []);
cfg.maskparameter = ft_getopt(cfg, 'maskparameter',[]);
cfg.linestyle     = ft_getopt(cfg, 'linestyle',    '-');
cfg.linewidth     = ft_getopt(cfg, 'linewidth',    0.5);
cfg.maskstyle     = ft_getopt(cfg, 'maskstyle',    'box');
cfg.channel       = ft_getopt(cfg, 'channel',      'all');
cfg.matrixside    = ft_getopt(cfg, 'matrixside',   'outflow');

ndata = numel(varargin);

% interactive plotting is not allowed with more than 1 input
% if ndata >1 && strcmp(cfg.interactive, 'yes')
%   error('interactive plotting is not supported with more than 1 input data set');
% end

%fixme rename matrixside and cohrefchannel in more meaningful options
if ischar(cfg.graphcolor)
    graphcolor = ['k' cfg.graphcolor];
elseif isnumeric(cfg.graphcolor)
    graphcolor = [0 0 0; cfg.graphcolor];
end

% check for linestyle being a cell-array, check it's length, and lengthen it if does not have enough styles in it
if ischar(cfg.linestyle)
    cfg.linestyle = {cfg.linestyle};
end

if ndata  > 1
    if (length(cfg.linestyle) < ndata ) && (length(cfg.linestyle) > 1)
        error('either specify cfg.linestyle as a cell-array with one cell for each dataset, or only specify one linestyle')
    elseif (length(cfg.linestyle) < ndata ) && (length(cfg.linestyle) == 1)
        tmpstyle = cfg.linestyle{1};
        cfg.linestyle = cell(ndata , 1);
        for idataset = 1:ndata
            cfg.linestyle{idataset} = tmpstyle;
        end
    end
end

% ensure that the input is correct, also backward compatibility with old data structures:
dtype = cell(ndata , 1);
for i=1:ndata
    varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'timelock', 'freq'});
    dtype{i}    = ft_datatype(varargin{i});
    
    % this is needed for correct treatment of graphcolor later on
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

if ndata >1,
    if ~all(strcmp(dtype{1}, dtype))
        error('input data are of different type; this is not supported');
    end
end
dtype  = dtype{1};
dimord = varargin{1}.dimord;
dimtok = tokenize(dimord, '_');


% set x/y/zparam defaults according to datatype and dimord
switch dtype
    case 'timelock'
        cfg.xparam = ft_getopt(cfg,  'xparam', 'time');
        cfg.yparam = ft_getopt(cfg,  'yparam', '');
        cfg.zparam = ft_getopt(cfg,  'zparam', 'avg');
    case 'freq'
        if sum(ismember(dimtok, 'time'))
            cfg.xparam = ft_getopt(cfg,  'xparam', 'time');
            cfg.yparam = ft_getopt(cfg,  'yparam', 'freq');%fixme
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
    for i=1:ndata
        varargin{i} = ft_timelockanalysis(tmpcfg, varargin{i});
    end
    dimord        = varargin{1}.dimord;
    dimtok        = tokenize(dimord, '_');
elseif strcmp(dtype, 'freq') && hasrpt,
    % this also deals with fourier-spectra in the input
    % or with multiple subjects in a frequency domain stat-structure
    % on the fly computation of coherence spectrum is not supported
    for i=1:ndata
        if isfield(varargin{i}, 'crsspctrm'),
            varargin{i} = rmfield(varargin{i}, 'crsspctrm');
        end
    end
    
    tmpcfg           = [];
    tmpcfg.trials    = cfg.trials;
    tmpcfg.jackknife = 'no';
    for i=1:ndata
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

% apply baseline correction
if ~strcmp(cfg.baseline, 'no')
    for i=1:ndata
        if strcmp(dtype, 'timelock') && strcmp(cfg.xparam, 'time')
            varargin{i} = ft_timelockbaseline(cfg, varargin{i});
        elseif strcmp(dtype, 'freq') && strcmp(cfg.xparam, 'time')
            varargin{i} = ft_freqbaseline(cfg, varargin{i});
        elseif strcmp(dtype, 'freq') && strcmp(cfg.xparam, 'freq')
            error('baseline correction is not supported for spectra without a time dimension');
        else
            warning('baseline correction not applied, please set cfg.xparam');
        end
    end
end

% handle the bivariate case

% check for bivariate metric with 'chan_chan' in the dimord
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;

% check for bivariate metric with a labelcmb
haslabelcmb = isfield(varargin{1}, 'labelcmb');

if (isfull || haslabelcmb) && isfield(varargin{1}, cfg.zparam)
    % a reference channel is required:
    if ~isfield(cfg, 'refchannel')
        error('no reference channel is specified');
    end
    
    % check for refchannel being part of selection
    if ~strcmp(cfg.refchannel,'gui')
        if (isfull      && ~any(ismember(varargin{1}.label, cfg.refchannel))) || ...
                (haslabelcmb && ~any(ismember(varargin{1}.labelcmb(:), cfg.refchannel)))
            error('cfg.refchannel is a not present in the (selected) channels)')
        end
    end
    
    % interactively select the reference channel
    if strcmp(cfg.refchannel, 'gui')
        error('cfg.refchannel = ''gui'' is not supported in ft_singleplotER');
    end
    
    for i=1:ndata
        if ~isfull,
            % convert 2-dimensional channel matrix to a single dimension:
            if isempty(cfg.matrixside)
                sel1 = strmatch(cfg.refchannel, varargin{i}.labelcmb(:,2), 'exact');
                sel2 = strmatch(cfg.refchannel, varargin{i}.labelcmb(:,1), 'exact');
            elseif strcmp(cfg.matrixside, 'outflow')
                sel1 = [];
                sel2 = strmatch(cfg.refchannel, varargin{i}.labelcmb(:,1), 'exact');
            elseif strcmp(cfg.matrixside, 'inflow')
                sel1 = strmatch(cfg.refchannel, varargin{i}.labelcmb(:,2), 'exact');
                sel2 = [];
            end
            fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.zparam);
            varargin{i}.(cfg.zparam) = varargin{i}.(cfg.zparam)([sel1;sel2],:,:);
            varargin{i}.label     = [varargin{i}.labelcmb(sel1,1);varargin{i}.labelcmb(sel2,2)];
            varargin{i}.labelcmb  = varargin{i}.labelcmb([sel1;sel2],:);
            varargin{i}           = rmfield(varargin{i}, 'labelcmb');
        else
            % general case
            sel               = match_str(varargin{i}.label, cfg.refchannel);
            siz               = [size(varargin{i}.(cfg.zparam)) 1];
            if strcmp(cfg.matrixside, 'inflow') || isempty(cfg.matrixside)
                %the interpretation of 'inflow' and 'outflow' depend on
                %the definition in the bivariate representation of the data  
                %data.(cfg.zparam) = reshape(mean(data.(cfg.zparam)(:,sel,:),2),[siz(1) 1 siz(3:end)]);
                sel1 = 1:siz(1);
                sel2 = sel;
                meandir = 2;
            elseif strcmp(cfg.matrixside, 'outflow')
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

% get physical min/max range of x
if strcmp(cfg.xlim,'maxmin')
    % find maxmin throughout all varargins:
    xmin = [];
    xmax = [];
    for i=1:ndata
        xmin = min([xmin varargin{i}.(cfg.xparam)]);
        xmax = max([xmax varargin{i}.(cfg.xparam)]);
    end
else
    xmin = cfg.xlim(1);
    xmax = cfg.xlim(2);
end

% get the index of the nearest bin
for i=1:ndata
    xidmin(i,1) = nearest(varargin{i}.(cfg.xparam), xmin);
    xidmax(i,1) = nearest(varargin{i}.(cfg.xparam), xmax);
end

%fixme do something with yparam here
%technically should not be defined for multiplotER, but can be defined (and
%use ft_selectdata to average across frequencies

cla
hold on;
colorlabels = [];

% plot each data set:
for i=1:ndata
    if isfield(varargin{1}, 'label')
        selchannel = ft_channelselection(cfg.channel, varargin{i}.label);
    elseif isfield(varargin{1}, 'labelcmb')
        selchannel = ft_channelselection(cfg.channel, unique(varargin{i}.labelcmb(:)));
    else
        error('the input data does not contain a label or labelcmb-field');
    end
    
    % make vector dat with one value for each channel
    dat    = varargin{i}.(cfg.zparam);
    xparam = varargin{i}.(cfg.xparam);
    
    % take subselection of channels
    % this works for bivariate data with labelcmb because at this point the
    % data has a label-field
    sellab = match_str(varargin{i}.label, selchannel);
    
    if ~isempty(cfg.yparam)
        if isfull
            dat = dat(sel1, sel2, ymin:ymax, xidmin(i):xidmax(i));
            dat = nanmean(nanmean(dat, meandir), 3);
            siz = size(dat);
            %fixmedat = reshape(dat, [siz(1:2) siz(4)]);
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
            siz(find(siz(1:2)==1)) = [];            
            dat = reshape(dat, siz);
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
        datmask = varargin{1}.(cfg.maskparameter)(sellab,:);
        datmask = datmask(xidmin(i):xidmax(i));
        maskdatavector = reshape(mean(datmask,1), [1 numel(xvector)]);
    else
        maskdatavector = [];
    end
    
    if ndata  > 1
        if ischar(graphcolor);        colorlabels = [colorlabels iname{i+1} '=' graphcolor(i+1) '\n'];
        elseif isnumeric(graphcolor); colorlabels = [colorlabels iname{i+1} '=' num2str(graphcolor(i+1,:)) '\n'];
        end
    end
    
    if ischar(graphcolor);        color = graphcolor(i+1);
    elseif isnumeric(graphcolor); color = graphcolor(i+1,:);
    end
    
    % update ymin and ymax for the current data set:
    if strcmp(cfg.ylim,'maxmin')
        if i==1
            ymin = [];
            ymax = [];
        end
        % select the channels in the data that match with the layout:
        ymin = min([ymin min(datavector)]);
        ymax = max([ymax max(datavector)]);
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

% set xlim and ylim:
xlim([xmin xmax]);
ylim([ymin ymax]);

% adjust mask box extents to ymin/ymax
ptchs = findobj(gcf,'type','patch');
for i = 1:length(ptchs)
    YData = get(ptchs(i),'YData');
    YData(YData == min(YData)) = ymin;
    YData(YData == max(YData)) = ymax;
    set(ptchs(i),'YData',YData);
end

if strcmp('yes',cfg.hotkeys)
    %  attach data and cfg to figure and attach a key listener to the figure
    set(gcf, 'keypressfcn', {@key_sub, xmin, xmax, ymin, ymax})
end

% make the figure interactive
if strcmp(cfg.interactive, 'yes')
    set(gcf, 'windowbuttonupfcn',     {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'windowbuttonupfcn'});
    set(gcf, 'windowbuttondownfcn',   {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'windowbuttondownfcn'});
    set(gcf, 'windowbuttonmotionfcn', {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg, varargin{:}}, 'event', 'windowbuttonmotionfcn'});
end

% create title text containing channel name(s) and channel number(s):
if length(sellab) == 1
    t = [char(cfg.channel) ' / ' num2str(sellab) ];
else
    t = sprintf('mean(%0s)', join(',', cfg.channel));
end
h = title(t,'fontsize', cfg.fontsize);

% set renderer if specified
if ~isempty(cfg.renderer)
    set(gcf, 'renderer', cfg.renderer)
end

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
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
% subfunction which is called by ft_select_channel in case cfg.refchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label, cfg, varargin)
cfg.refchannel = label;
fprintf('selected cfg.refchannel = ''%s''\n', cfg.refchannel);
p = get(gcf, 'position');
f = figure;
set(f, 'position', p);
ft_singleplotER(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotER(range, cfg, varargin)
cfg.comment = 'auto';
cfg.yparam = [];
cfg.xlim = range(1:2);
if isfield(cfg, 'showlabels')
  % this is not allowed in topoplotER
  cfg = rmfield(cfg, 'showlabels');
end
fprintf('selected cfg.xlim = [%f %f]\n', cfg.xlim(1), cfg.xlim(2));
p = get(gcf, 'position');
f = figure;
set(f, 'position', p);
ft_topoplotER(cfg, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata, varargin)
xlimits = xlim;
ylimits = ylim;
incr_x = abs(xlimits(2) - xlimits(1)) /10;
incr_y = abs(ylimits(2) - ylimits(1)) /10;

% TRANSLATE by 10%
if length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:},'control') && strcmp(eventdata.Key,'leftarrow') 
    xlim([xlimits(1)+incr_x xlimits(2)+incr_x])
elseif length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:},'control') && strcmp(eventdata.Key,'rightarrow')  
    xlim([xlimits(1)-incr_x xlimits(2)-incr_x])
elseif length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:},'control') && strcmp(eventdata.Key,'uparrow')
    ylim([ylimits(1)-incr_y ylimits(2)-incr_y])
elseif length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:},'control') && strcmp(eventdata.Key,'downarrow')
   ylim([ylimits(1)+incr_y ylimits(2)+incr_y])
% ZOOM by 10%
elseif strcmp(eventdata.Key,'leftarrow') 
    xlim([xlimits(1)-incr_x xlimits(2)+incr_x])
elseif strcmp(eventdata.Key,'rightarrow') 
    xlim([xlimits(1)+incr_x xlimits(2)-incr_x])
elseif strcmp(eventdata.Key,'uparrow')
    ylim([ylimits(1)-incr_y ylimits(2)+incr_y])
elseif strcmp(eventdata.Key,'downarrow') 
    ylim([ylimits(1)+incr_y ylimits(2)-incr_y])
% resort to minmax of data for x-axis and y-axis
elseif strcmp(eventdata.Key,'m')
    xlim([varargin{1} varargin{2}])
    ylim([varargin{3} varargin{4}])
end


