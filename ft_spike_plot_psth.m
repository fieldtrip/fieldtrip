function [hdl] = ft_spike_plot_psth(cfg,psth)

% FT_SPIKE_PLOT_PSTH makes a bar plot of PSTH structure (output from FT_SPIKE_PSTH) 
% with error bars.
%
%	Inputs:
%		PSTH typically is a structure from FT_SPIKE_PSTH.
%
%   Configurations (cfg):
%
%   cfg.latency          = [begin end] in seconds, 'maxperiod' (default), 'prestim'(t<=0), or
%                          'poststim' (t>=0).     
%   cfg.errorbars        = 'no', 'std', 'sem' (default), 'conf95%' (requires statistic toolbox, 
%                          according to student-T distribution), 'var'
%   cfg.spikechannel     = string or index of single spike channel to trigger on (default = 1)
%                          Only one spikechannel can be plotted at a time.
%   cfg.ylim             = [min max] or 'auto' (default)
%                          If 'standard', we plot from 0 to 110% of maximum plotted value);
%   Outputs:
%		HDL.avg 		 = figure handle for the bar plot, psth average. Use SET and GET to access.
%		HDL.var 		 = figure handle for the error lines. Use GET and SET to access.
%

% Martin Vinck (C) 2010.

defaults.latency      = {'maxperiod'}; 
defaults.errorbars    = {'sem' 'std' 'conf95%', 'no' 'var'}; 
defaults.spikechannel = {1}; 
defaults.ylim         = {'auto'};

cfg = ft_spike_sub_defaultcfg(cfg,defaults);% change prefix

% select the latencies
if ischar(cfg.latency)
  if strcmp(cfg.latency, 'maxperiod')
      cfg.latency = [min(psth.time) max(psth.time)];
  elseif strcmp(cfg.latency, 'prestim')
      cfg.latency = [min(psth.time) 0];
  elseif strcmp(cfg.latency, 'poststim')    
      cfg.latency = [0 max(psth.time)];
  end
elseif ~isnumeric(cfg.latency)||length(cfg.latency)~=2
  error('MATLAB:ft_spike_plot_psth:cfg:latency',...
    'cfg.latency should be "maxperiod", "prestim", "poststim" or 1-by-2 numerical vector');
end

if cfg.latency(1)>=cfg.latency(2), error('MATLAB:spikestation:plot_psth:incorrectLatencyWindow',...
    'cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end

% check whether the time window fits with the data
if cfg.latency(1) < min(psth.time), cfg.latency(1) = min(psth.time); 
  warning('MATLAB:ft_spike_plot_psth:incorrectLatencyWindow',...
          'Correcting begin latency of averaging window');
end
if (cfg.latency(2) > max(psth.time)), cfg.latency(2) = max(psth.time);
  warning('MATLAB:ft_spike_plot_psth:incorrectLatencyWindow',...
          'Correcting begin latency of averaging window');
end

% get the spikechannels
cfg.channel = ft_channelselection(cfg.spikechannel, psth.label);
spikesel    = match_str(psth.label, cfg.channel);
nUnits      = length(spikesel); 
if nUnits~=1, error('MATLAB:ft_spike_plot_psth:cfg:spikechannel:wrongInput',...
      'You selected more or less than one spikechannel by means of cfg.spikechannel'); 
end

% select the timepoints within the latencies
timeSel = psth.time>=cfg.latency(1) & psth.time <= cfg.latency(2);
if isempty(timeSel), error('MATLAB:ft_spike_plot_psth:cfg:latency',...
    'no time points selected, please change cfg.latency or inspect input PSTH');
end

% plot the average psth
psthHdl = bar(psth.time(timeSel),psth.avg(spikesel,timeSel),'k');
set(psthHdl,'BarWidth', 1)

% check if fields .dof and .var are correct given cfg.errorbars
if ~strcmp(cfg.errorbars,'no')
  if ~isfield(psth,'var') || ~any(isfinite(psth.var(spikesel,:))) 
     error('MATLAB:ft_spike_plot_psth:cfg:var',...
    'PSTH should contain field .var that contain numbers. If you do not want the variance',...
    'please specify cfg.errorbars = "no"');
  end
  if ~ (strcmp(cfg.errorbars,'std') || strcmp(cfg.errorbars,'var'))
    if ~isfield(psth,'dof') || ~any(isfinite(psth.dof(spikesel,:)))
       error('MATLAB:ft_spike_plot_psth:cfg:dof',...
      'psth should contain field dof that contains numbers. Use cfg.errorbars = "no" if you do not',...
      'want the variance plotted');
    end
  end
end

% add standard error bars to it
x = [psth.time(timeSel);psth.time(timeSel)]; % create x doublets
if strcmp(cfg.errorbars, 'sem')
    err = sqrt(psth.var(spikesel,timeSel)./psth.dof(timeSel));
elseif strcmp(cfg.errorbars, 'std')
    err = sqrt(psth.var(spikesel,timeSel));
elseif strcmp(cfg.errorbars, 'var')
    err = psth.var(spikesel,timeSel);
elseif strcmp(cfg.errorbars, 'conf95%')
    % use a try statement just in case the statistics toolbox doesn't work.
    try
        tCrit = tinv(0.975,psth.dof(timeSel));
        err = tCrit.*sqrt(psth.var(spikesel,timeSel)./psth.dof(timeSel)); % assuming normal distribution, SHOULD BE REPLACED BY STUDENTS-T!
    end
else
    err = 0;
end
y = psth.avg(spikesel,timeSel) + err; % create the error values

% plot the errorbars as line plot on top of the bar plot
if ~strcmp(cfg.errorbars,'no')
  hold on
  y = [zeros(1,length(y));y]; % let lines run from 0 to error values
  psthSemHdl = plot(x,y,'k'); % plot the error bars  
end

% create labels
xlabel('bin centers (sec)')
try
  ylabel(psth.cfg.outputunit)
catch
  ylabel('firing activity')
end

% set the axis
if strcmp(cfg.ylim, 'auto'), cfg.ylim = [0 max(y(:))*1.1+eps]; end
if ~isnumeric(cfg.ylim)||length(cfg.ylim)~=2 || cfg.ylim(2)<=cfg.ylim(1)
    error('MATLAB:spikestation:plot_psth:cfg:ylim',...
    'cfg.ylim should be "auto" or ascending order 1-by-2 vector in seconds');
end
set(gca,'YLim', cfg.ylim, 'XLim', cfg.latency)
set(gca,'TickDir','out', 'Box', 'off')

% store the handles
hdl.avg    = psthHdl;
if ~strcmp(cfg.errorbars,'no'), hdl.error = psthSemHdl;  end

% as usual, make sure that panning and zooming does not distort the y limits
set(zoom,'ActionPostCallback',{@mypostcallback,cfg.ylim,cfg.latency});
set(pan,'ActionPostCallback',{@mypostcallback,cfg.ylim,cfg.latency});

function [] = mypostcallback(fig,evd,limY,limX)

% get the x limits and reset them
ax = evd.Axes;
xlim = get(ax(1), 'XLim');
ylim = get(ax(1), 'YLim');
if limX(1)>xlim(1), xlim(1) = limX(1); end
if limX(2)<xlim(2), xlim(2) = limX(2); end      
if limY(1)>ylim(1), ylim(1) = limY(1); end
if limY(2)<ylim(2), ylim(2) = limY(2); end      

set(ax(1), 'XLim',xlim)
set(ax(1), 'YLim',ylim);  


