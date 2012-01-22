function [cfg] = ft_spike_plot_raster(cfg, spike)

% FT_SPIKE_PLOT_RASTER makes a raster plot of spike-trains and allows for a
% spike-density or a PSTH plot on top.
%
% Use as
%   ft_spike_plot_raster(cfg, spike)
% 
% The input spike data structure should be organised as the spike or the raw
% datatype, obtained from FT_SPIKE_MAKETRIALS or FT_PREPROCESSING (in that case
% conversion is done within the function).
%
% Configuration options related to selection of spike channel and trials and latencies
%   cfg.spikechannel     =  see FT_CHANNELSELECTION for details
%   cfg.latency          =  [begin end]` in seconds, 'maxperiod' (default), 'minperiod',
%                           'prestim' (all t<=0), or 'poststim' (all t>=0).
%   cfg.linewidth        =  number indicating the width of the lines (default = 1);
%   cfg.cmapneurons      =  'auto' (default), or nUnits-by-3 matrix, or cell array with
%                           color strings (e.g., {'k' 'r' 'b'}). If 'auto', we are using a
%                           private colormap that has good visibility on white background.
%   cfg.spikelength      =  number >0 and <=1 indicating the length of the spike. If
%                           cfg.spikelength = 1, then no space will be left between
%                           subsequent rows representing trials (row-unit is 1).
%
% Configuration options related to additionally plotting the TOPDATA
%   cfg.topdata          =  output structure from FT_SPIKE_PSTH
%                           or FT_SPIKEDENSITY or FT_TIMELOCKANALYSIS. See those functions for more
%                           info.
%   cfg.topplotsize      =  number ranging from 0 to 1, indicating the proportion of the
%                           rasterplot that the top plot will take (e.g., with 0.7 the top
%                           plot will be 70% of the rasterplot in size). Default = 0.5.
%   cfg.topplotfunc      =  'bar' (default) or 'line'.
%
% Outputs:
%   hdl.top              =  handle of the plot of the topdata (psth or sdf)
%   hdl.raster           =  handle for each individual spike.
%   hdl.ax               =  axes handles: ax(1) the rasterplot and ax(2) the topplot.
%   hdl.cfg

% Copyright (C) 2010, Martin Vinck; F.C. Donders Centre for Neuroimaging; University of Amsterdam
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% check if input spike structure is indeed a spike structure
spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

% get the default options
cfg.spikechannel = ft_getopt(cfg, 'spikechannel', 'all');
cfg.trials       = ft_getopt(cfg, 'trials', 'all');
cfg.latency      = ft_getopt(cfg,'latency','maxperiod');
cfg.linewidth    = ft_getopt(cfg,'linewidth', 1);
cfg.cmapneurons  = ft_getopt(cfg,'cmapneurons', 'auto');
cfg.spikelength  = ft_getopt(cfg,'spikelength', 0.9);
cfg.topdata      = ft_getopt(cfg,'topdata', []);
cfg.topplotsize  = ft_getopt(cfg,'topplotsize', 0.5);
cfg.topplotfunc  = ft_getopt(cfg,'topplotfunc', 'bar');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'latency', {'char', 'doublevector'});
cfg = ft_checkopt(cfg,'trials', {'char', 'doublevector', 'logical'}); 
cfg = ft_checkopt(cfg,'linewidth', 'doublescalar');
cfg = ft_checkopt(cfg,'cmapneurons', {'char', 'double', 'cell'});
cfg = ft_checkopt(cfg,'spikelength', 'doublescalar');
cfg = ft_checkopt(cfg,'topdata', {'struct', 'empty'});
cfg = ft_checkopt(cfg,'topplotsize', 'doublescalar');
cfg = ft_checkopt(cfg,'topplotfunc', 'char', {'bar', 'line'});

% check which features should be present in the rasterplot and psth
if ~isempty(cfg.topdata)
  doTopData = true;
  topData = cfg.topdata;
  cfg = rmfield(cfg, 'topdata');
  topData = ft_checkdata(topData,'datatype', 'timelock', 'feedback', 'yes');
else
  doTopData = false;
end

% get the spikechannels
cfg.spikechannel = ft_channelselection(cfg.spikechannel, spike.label);
spikesel    = match_str(spike.label, cfg.spikechannel);
nUnits   = length(spikesel); % number of spike channels
if nUnits==0, error('MATLAB:ft_spike_plot_raster:cfg:spikechannel:noSpikeChanSelected',...
    'No spikechannel selected by means of cfg.spikechannel');
end

% get the number of trials or change DATA according to cfg.trials
nTrialsOrig = size(spike.trialtime,1);
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:nTrialsOrig;
elseif islogical(cfg.trials)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrialsOrig, error('MATLAB:ft_spike_plot_raster:cfg:trials:maxExceeded',...
    'maximum trial number in cfg.trials should not exceed length of DATA.trial')
end
if isempty(cfg.trials),
  errors('MATLAB:ft_spike_plot_raster:cfg:trials','No trials were selected in cfg.trials');
end
nTrials = length(cfg.trials);

% determine the duration of each trial
begTrialLatency = spike.trialtime(cfg.trials,1);
endTrialLatency = spike.trialtime(cfg.trials,2);

% select the latencies
if strcmp(cfg.latency,'minperiod')
  cfg.latency = [max(begTrialLatency) min(endTrialLatency)];
elseif strcmp(cfg.latency,'maxperiod')
  cfg.latency = [min(begTrialLatency) max(endTrialLatency)];
elseif strcmp(cfg.latency,'prestim')
  cfg.latency = [min(begTrialLatency) 0];
elseif strcmp(cfg.latency,'poststim')
  cfg.latency = [0 max(endTrialLatency)];
elseif ~isrealvec(cfg.latency)||length(cfg.latency)~=2
  error('MATLAB:fieldtrip:spikerate:cfg:latency',...
    'cfg.latency should be "max", "min", "prestim", "poststim" or 1-by-2 numerical vector');
end
if cfg.latency(1)>cfg.latency(2),
  error('MATLAB:ft_spike_plot_raster:cfg:latency:wrongOrder',...
    'cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end
% check whether the time window fits with the data
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency);
  warning('MATLAB:ft_spike_plot_raster:cfg:latency:correctBeg',...
    'Correcting begin latency of averaging window');
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('MATLAB:ft_spike_plot_raster:cfg:latency:correctEnd',...
    'Correcting begin latency of averaging window');
end
trialDur = cfg.latency(2) - cfg.latency(1);

% delete trials that are not in our window
overlaps      = endTrialLatency>=cfg.latency(1) & begTrialLatency<=cfg.latency(2);
trialSel      = overlaps(:);
cfg.trials    = cfg.trials(trialSel); %update the trial selection
if isempty(cfg.trials),
  warning('MATLAB:ft_spike_plot_raster:cfg:trials','No trials were selected');
end
nTrials         = length(cfg.trials);

% create the data that should be plotted
[unitX,unitY] = deal(cell(1,nUnits));

% find the spiketimes per neuron and make one vector of them with y value per element
for iUnit = 1:nUnits
  unitIndx       = spikesel(iUnit);
  latencySel     = spike.time{unitIndx}>=cfg.latency(1) & spike.time{unitIndx} <= cfg.latency(2);
  isInTrials     = ismember(spike.trial{unitIndx},cfg.trials);
  unitX{iUnit}   = spike.time{unitIndx}(isInTrials(:) & latencySel(:));
  unitY{iUnit}   = spike.trial{unitIndx}(isInTrials(:) & latencySel(:));
end

if ~isscalar(cfg.spikelength) || cfg.spikelength<=0 || cfg.spikelength>1
  error('MATLAB:ft_spike_plot_raster:cfg:spikelength:unknownOption',...
    'cfg.spikelength should be a single number >0 and <=1. 1 row (1 trial) = 1');
end

if ~isscalar(cfg.topplotsize) || cfg.topplotsize<=0 || cfg.topplotsize>1
  error('MATLAB:ft_spike_plot_raster:cfg:topplotsize:unknownOption',...
    'cfg.topplotsize should be a single number >0 and <=1. 0.7 = 70%');
end

% start plotting the rasterplot
for iUnit = 1:nUnits
  
  % create the x and y value to be plotted
  x = [unitX{iUnit}(:)'; unitX{iUnit}(:)']; % create x duplicates, spike times
  y = unitY{iUnit}(:); % these are the trial values
  y = [y' - cfg.spikelength/2; y' + cfg.spikelength/2];
  
  % process the color information for the different units
  if strcmp(cfg.cmapneurons, 'auto'), cfg.cmapneurons = colormap_cgbprb(nUnits); end
  if isrealmat(cfg.cmapneurons) && all(size(cfg.cmapneurons) ./ [nUnits 3])
    color = cfg.cmapneurons(iUnit,:);
  elseif iscell(cfg.cmapneurons) && length(cfg.cmapneurons)==nUnits
    color = cfg.cmapneurons{iUnit};
  else
    error('MATLAB:ft_spike_plot_raster:cfg:cmapneurons:unknownOption',...
      'cfg.cmapneurons should be nUnits-by-3 matrix or 1-by-nUnits cell or "auto"');
  end
  
  % create axes for the rasterplot, all go to the same position, so do this for unit 1
  if iUnit==1
    ax(1) = newplot;
    posOrig = get(ax(1), 'Position'); % original position for a standard plot
    pos     = posOrig;
    if doTopData % if topdata, we leave space on top
      posRaster   = [pos(1:3) pos(4)*(1-cfg.topplotsize)]; % and decrease the size of width and height
    else
      posRaster   = pos;
    end
    set(ax(1), 'ActivePositionProperty', 'position', 'Position', posRaster)
  end
  
  % make the raster plot and hold on for the next plots
  rasterHdl = plot(x, y,'linewidth', cfg.linewidth,'Color', color);
  set(ax(1),'NextPlot', 'add')
  set(ax(1),'Box', 'off')
end

% create the labels for the first axes
xlabel('time (sec)')
ylabel('Trial Number')
axis ij

% plot the top data
if doTopData
  
  % match on the basis of labels, and specify this in the documentary
  unitIndx = find(ismember(topData.label,spike.label(spikesel)));
  if sum(unitIndx)<nUnits, error('MATLAB:ft_spike_plot_raster:topData:missingLabel',...
      'topData.label should contain all label of selected units');
  end
  
  % select timepoints  and get the data to be plotted
  binSel = topData.time>=cfg.latency(1) & topData.time<=cfg.latency(2);
  y      = topData.avg(unitIndx,binSel);
  
  % create the position for the topplot and axes with this position
  posTopPlot = [0 pos(4)*(1-cfg.topplotsize)  0 0] + pos.*[1 1 1 cfg.topplotsize];
  ax(2)  = axes('Units', get(ax(1), 'Units'), 'Position', posTopPlot,...
    'Parent', get(ax(1), 'Parent'));
  
  % make the bar or the line plot
  if strcmp(cfg.topplotfunc,'bar')
    avgHdl  = feval(cfg.topplotfunc,topData.time(binSel),y', 'Stack');
    for iUnit = 1:nUnits
      set(avgHdl(iUnit),'FaceColor',cfg.cmapneurons(iUnit,:),'EdgeColor', cfg.cmapneurons(iUnit,:));
    end
    if strcmp(cfg.topplotfunc,'bar'),set(avgHdl,'BarWidth', 1); end
  else
    x = topData.time(binSel);
    x = repmat(x(:),1,nUnits); % it puts multiple neurons here? check>
    avgHdl  = feval(cfg.topplotfunc,x,y');
    for iUnit = 1:nUnits
      set(avgHdl(iUnit),'Color',cfg.cmapneurons(iUnit,:));
    end
  end
  
  % modify the axes
  set(ax(2),'YAxisLocation', 'right') % swap y axis location
  set(ax(2),'XTickLabel', {}) % remove ticks and labels for x axis
  set(ax(2),'XTick', [])
  set(ax(2), 'Box', 'off') % turn the box off
  
  % change the axis settings
  try
    ylabel(topData.cfg.outputunit)
  catch
    ylabel('Firing Rate')
  end
end

% set the limits for the axis
set(ax,'XLim', [cfg.latency])
set(ax(1), 'YLim', [0.5 nTrialsOrig+0.5]); % number of trials
set(ax,'TickDir','out') % put the tickmarks outside

% collect the handles
hdl.raster = rasterHdl;
if doTopData,                 hdl.top = avgHdl;       end
hdl.ax     = ax;

% now link the axes, constrain zooming and keep ticks intact
limX       = [cfg.latency];
limY       = get(ax,'YLim');
if ~iscell(limY), limY = {limY}; end

% constrain the zooming and zoom psth together with the jpsth, remove ticklabels jpsth
set(zoom,'ActionPostCallback',{@mypostcallback,ax,limX,limY});
set(pan,'ActionPostCallback',{@mypostcallback,ax,limX,limY});
hdl.cfg = cfg;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous spike

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = mypostcallback(fig,evd,ax,lim,ylim)
currentAxes = evd.Axes;
indx = find(currentAxes==ax);

% keep the y axis within boundaries
origLim = ylim{indx};
ylim = get(ax(indx), 'YLim');
if origLim(1)>ylim(1), ylim(1) = origLim(1); end
if origLim(2)<ylim(2), ylim(2) = origLim(2); end
set(ax(indx), 'YLim', ylim)

% keep the x axis within boundaries and reset for both
xlim = get(ax(indx), 'XLim');
if lim(1)>xlim(1), xlim(1) = lim(1); end
if lim(2)<xlim(2), xlim(2) = lim(2); end
set(ax,'XLim', xlim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y] = polygon(x,mn,sm,multiplier)
if nargin < 4
  multiplier = 1;
end
X = [x x(end:-1:1) x(1)];
up = mn + multiplier*sm;
down = mn - multiplier*sm;
Y = [down up(end:-1:1) down(1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [colorspace] = colormap_cgbprb(N)
% COLORMAP_SEPARATION returns a color map with well separated colors that does not include
% white, yellow or orange since these are not well visible in white plots.
%
% Inputs: N is the number of colors desired
% Outputs: COLORSPACE is an N-by-3 colorspace.
% Chosen colors are green, cyan, blue, purple, red and black
% Script was adapted from varycolors.m from FEX
% Copyright (c) 2009, Martin Vinck

if nargin<1
  error('specify the number of colors that you want')
end

% create a well visible color space
%        green     cyan  blue    purple    red     black
colors = [[0 1 0]; [0 1 1] ; [0 0 1] ; [1 0 1]; [1 0 0]; [0 0 0]];
order  = [1 5 3 2 4 6]; % our preference order when plotting different lines
rm{2} = [2:5];
rm{3} = [2 4 5];
rm{4} = [2 4];
rm{5} = [4];

if N==1
  colorspace = colors(end,:);
elseif N>1&&N<6
  colors(rm{N},:) = [];
  colorspace = colors;
else
  n      = floor(N/5)*ones(1,5);
  rest   = mod(N,5);
  order  = [1 5 3 2 4];
  % if we have some rest, add them starting oppositly
  n(order(1:rest)) = n(order(1:rest)) + 1;
  
  colorspace = [];
  for iColor = 1:5
    corStart  = colors(iColor,:);
    corEnd    = colors(iColor+1,:);
    dim       = corStart~=corEnd;
    subColors = corStart(ones(n(iColor)+1,1),:);
    if iColor>1
      subColors(:,dim) = linspace(corStart(dim),corEnd(dim),n(iColor)+1);
      subColors(1,:) = [];
    else
      subColors(1:end-1,dim) =   linspace(corStart(dim),corEnd(dim),n(iColor));
      subColors(end,:) = [];
    end
    colorspace = [colorspace;subColors];
  end
end
