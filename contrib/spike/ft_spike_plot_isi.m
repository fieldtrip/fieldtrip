function [cfg] = ft_spike_plot_isi(cfg, isih)

% FT_SPIKE_PLOT_ISI makes an inter-spike-interval bar plot.
%
% Use as
%   ft_spike_plot_isi(cfg, isih)
%
% Inputs:
%   ISIH is the output from FT_SPIKE_ISIHIST
%
% Configurations:
%   cfg.spikechannel     = string or index or logical array to to select 1 spike channel.
%                          (default = 1).
%   cfg.ylim             = [min max] or 'auto' (default)
%                          If 'auto', we plot from 0 to 110% of maximum plotted value);
%   cfg.plotfit          = 'yes' (default) or 'no'. This requires that when calling
%                          FT_SPIKESTATION_ISI, cfg.gammafit = 'yes'.
%
% Outputs:
%   hdl.fit              = handle for line fit. Use SET and GET to access.
%   hdl.isih             = handle for bar isi histogram. Use SET and GET to access.

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
ft_preamble provenance isih


% get the default options
cfg.spikechannel = ft_getopt(cfg,'spikechannel', 'all');
cfg.ylim         = ft_getopt(cfg,'ylim', 'auto');
cfg.plotfit      = ft_getopt(cfg,'plotfit', 'no');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'ylim', {'char','ascendingdoublebivector'});
cfg = ft_checkopt(cfg,'plotfit', 'char', {'yes', 'no'});

% get the spikechannels
cfg.channel = ft_channelselection(cfg.spikechannel, isih.label);
spikesel    = match_str(isih.label, cfg.channel);
nUnits      = length(spikesel); % number of spike channels
if nUnits~=1, error('One spikechannel should be selected by means of cfg.spikechannel'); end

% plot the average isih
isiHdl = bar(isih.time,isih.avg(spikesel,:),'k');
set(isiHdl,'BarWidth', 1)

if strcmp(cfg.plotfit,'yes')
  
  if ~isfield(isih,'gamfit'), error('isih.gamfit should be present when cfg.plotfit = yes'); end
  
  % generate the probabilities according to the gamma model
  pGamma = gampdf(isih.time, isih.gamfit(spikesel,1), isih.gamfit(spikesel,2));
  sel    = isfinite(pGamma);
  pGamma = pGamma(sel)./sum(pGamma(sel));
  
  % scale the fit in the same way (total area is equal)
  sumIsih = sum(isih.avg(spikesel,:),2);
  pGamma  = pGamma*sumIsih;
  
  % plot the fit
  hold on
  fitHdl = plot(isih.time(sel), pGamma,'r');
else
  pGamma = 0;
end

try
  if strcmp(isih.cfg.outputunit, 'proportion')
    ylabel('probability')
  elseif strcmp(isih.cfg.outputunit, 'spikecount')
    ylabel('spikecount')
  end
catch
  ylabel('intensity');
end
    
xlabel('bin edges (sec)')

% set the x axis
set(gca,'XLim', [min(isih.time) max(isih.time)])

% set the y axis
if strcmp(cfg.ylim, 'auto')
  y = [isih.avg(spikesel,:) pGamma];
  cfg.ylim = [0 max(y(:))*1.1+eps];
end
set(gca,'YLim', cfg.ylim)
set(gca,'TickDir','out')

% store the handles
if strcmp(cfg.plotfit,'yes'),     H.fit    = fitHdl;    end
H.isih   = isiHdl;

% as usual, make sure that panning and zooming does not distort the y limits
set(zoom,'ActionPostCallback',{@mypostcallback,cfg.ylim,[min(isih.time) max(isih.time)]});
set(pan,'ActionPostCallback',{@mypostcallback,cfg.ylim,[min(isih.time) max(isih.time)]});

% do the general cleanup and bookkeeping at the end of the function

ft_postamble previous isih
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

