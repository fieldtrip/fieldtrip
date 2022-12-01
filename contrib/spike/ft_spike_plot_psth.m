function [cfg] = ft_spike_plot_psth(cfg, psth)

% FT_SPIKE_PLOT_PSTH makes a bar plot of PSTH structure with error bars.
%
% Use as
%   ft_spike_plot_psth(cfg, psth)
% 
%	Inputs:
%		PSTH typically is a structure from FT_SPIKE_PSTH.
%
% Configurations:
%   cfg.latency          = [begin end] in seconds, 'maxperiod' (default), 'prestim'(t<=0), or
%                          'poststim' (t>=0).
%   cfg.errorbars        = 'no', 'std', 'sem' (default), 'conf95%' (requires statistic toolbox,
%                          according to student-T distribution), 'var'
%   cfg.spikechannel     = string or index of single spike channel to trigger on (default = 1)
%                          Only one spikechannel can be plotted at a time.
%   cfg.ylim             = [min max] or 'auto' (default)
%                          If 'standard', we plot from 0 to 110% of maximum plotted value);
% Outputs:
%	  cfg.hdl.avg              = figure handle for the bar plot, psth average.
%	  cfg.hdl.var              = figure handle for the error lines. 
%
% See also FT_SPIKE_PSTH

% Copyright (C) 2010, Martin Vinck
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
ft_preamble provenance psth


% check the psth datatype
psth = ft_checkdata(psth, 'datatype', 'timelock', 'feedback', 'yes');

% get the default options
cfg.spikechannel = ft_getopt(cfg, 'spikechannel', 'all');
cfg.latency      = ft_getopt(cfg,'latency','maxperiod');
cfg.errorbars    = ft_getopt(cfg,'errorbars', 'sem');
cfg.ylim         = ft_getopt(cfg,'ylim', 'auto');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'latency', {'char', 'ascendingdoublebivector'});
cfg = ft_checkopt(cfg,'errorbars', 'char', {'sem', 'std', 'conf95%', 'no', 'var'});
cfg = ft_checkopt(cfg,'ylim', {'char', 'ascendingdoublebivector'});

cfg = ft_checkconfig(cfg, 'allowed', {'spikechannel', 'latency', 'errorbars', 'ylim'});

% select the latencies
if ischar(cfg.latency)
  if strcmp(cfg.latency, 'maxperiod')
    cfg.latency = [min(psth.time) max(psth.time)];
  elseif strcmp(cfg.latency, 'prestim')
    cfg.latency = [min(psth.time) 0];
  elseif strcmp(cfg.latency, 'poststim')
    cfg.latency = [0 max(psth.time)];
  end
end

% check whether the time window fits with the data
if cfg.latency(1) < min(psth.time), cfg.latency(1) = min(psth.time);
  warning('Correcting begin latency of averaging window');
end
if (cfg.latency(2) > max(psth.time)), cfg.latency(2) = max(psth.time);
  warning('Correcting begin latency of averaging window');
end

% get the spikechannels
cfg.spikechannel = ft_channelselection(cfg.spikechannel, psth.label);
spikesel    = match_str(psth.label, cfg.spikechannel);
nUnits      = length(spikesel);
if nUnits~=1, error('You selected more or less than one spikechannel by means of cfg.spikechannel'); end

% select the timepoints within the latencies
timeSel = psth.time>=cfg.latency(1) & psth.time <= cfg.latency(2);
if isempty(timeSel), error('no time points selected, please change cfg.latency or inspect input PSTH');end

% plot the average psth
psthHdl = bar(psth.time(timeSel),psth.avg(spikesel,timeSel),'k');
set(psthHdl,'BarWidth', 1)

% check if fields .dof and .var are correct given cfg.errorbars
if ~strcmp(cfg.errorbars,'no')
  if ~isfield(psth,'var') || ~any(isfinite(psth.var(spikesel,:)))
    error('PSTH should contain field .var with numbers');
  end
  if ~ (strcmp(cfg.errorbars,'std') || strcmp(cfg.errorbars,'var'))
    if ~isfield(psth,'dof') || ~any(isfinite(psth.dof(spikesel,:)))
      error('psth should contain field dof that contains numbers.'); 
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
  catch
    err = 0;
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
  ylabel('firing intensity')
end

% set the axis
if strcmp(cfg.ylim, 'auto'), cfg.ylim = [0 max(y(:))*1.1+eps]; end
if cfg.ylim(2)<max(y(:))
  warning('correcting cfg.ylim(2) as it clips the data'); 
  cfg.ylim(2) = max(y(:))*1.1+eps; 
end

set(gca,'YLim', cfg.ylim, 'XLim', cfg.latency)
set(gca,'TickDir','out', 'Box', 'off')

% store the handles
cfg.hdl.avg    = psthHdl;
if ~strcmp(cfg.errorbars,'no'), cfg.hdl.var = psthSemHdl;  end

% as usual, make sure that panning and zooming does not distort the y limits
set(zoom,'ActionPostCallback',{@mypostcallback,cfg.ylim,cfg.latency});
set(pan,'ActionPostCallback',{@mypostcallback,cfg.ylim,cfg.latency});

% do the general cleanup and bookkeeping at the end of the function

ft_postamble previous psth
ft_postamble provenance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

