function [cfg] = ft_spike_plot_jpsth(cfg, jpsth)

% FT_SPIKE_PLOT_JPSTH makes a plot from JPSTH structure.
%
% Use as
%   ft_spike_plot_jpsth(cfg, jpsth)
%
% Inputs:
%   JPSTH must be the output structure from FT_SPIKE_JPSTH and contain the
%   field JPSTH.avg. If cfg.psth = 'yes', the field JPSTH.psth must be
%   present as well.
%
% General configurations:
%   cfg.channelcmb  = string or index of single channel combination to trigger on.
%       See SPIKESTATION_FT_SUB_CHANNELCOMBINATION for details.
%   cfg.psth        = 'yes' (default) or 'no'. Plot PSTH with JPSTH if 'yes';
%   cfg.latency     = [begin end] in seconds or 'max' (default), 'prestim' or 'poststim';
%   cfg.colorbar    = 'yes' (default) or 'no'
%   cfg.colormap    =  N-by-3 colormap (see COLORMAP). or 'auto' (default,hot(256))
%   cfg.interpolate      = integer (default = 1), we perform interpolating
%                          with extra number of spacings determined by
%                          cfg.interpolate. For example cfg.interpolate = 5
%                          means 5 times more dense axis.
%   cfg.window      = 'string' or N-by-N matrix
%     'no':           apply no smoothing
%     ' gausswin'     use a Gaussian smooth function
%     ' boxcar'       use a box-car to smooth
%   cfg.gaussvar    =  variance  (default = 1/16 of window length in sec).
%   cfg.winlen      =  window length in seconds (default = 5*binwidth).
%     length of our window is 2*round*(cfg.winlen/binwidth)
%     where binwidth is the binwidth of the jpsth (jpsth.time(2)-jpsth.time(1)).
%
% See also FT_SPIKE_JPSTH, FT_SPIKE_PSTH

% FIXME: extend the windowing functions a bit

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
ft_preamble provenance jpsth


% get the default options
cfg.channelcmb  = ft_getopt(cfg,'channelcmb', 'all');
cfg.psth        = ft_getopt(cfg,'psth', 'yes');
cfg.latency     = ft_getopt(cfg,'latency','maxperiod');
cfg.colorbar    = ft_getopt(cfg,'colorbar', 'yes');
cfg.colormap    = ft_getopt(cfg,'colormap', jet(256));
cfg.interpolate = ft_getopt(cfg,'interpolate', 1);
cfg.window      = ft_getopt(cfg,'window', 'no');
cfg.winlen      = ft_getopt(cfg,'winlen', 5*(mean(diff(jpsth.time))));
cfg.gaussvar    = ft_getopt(cfg,'gaussvar', (cfg.winlen/4).^2);

% ensure that the options are valid
cfg = ft_checkopt(cfg,'channelcmb', {'char', 'cell'});
cfg = ft_checkopt(cfg,'psth', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'latency', {'char', 'ascendingdoublebivector'});
cfg = ft_checkopt(cfg,'colorbar', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'colormap','double');
cfg = ft_checkopt(cfg,'interpolate', 'doublescalar');
cfg = ft_checkopt(cfg,'window','char', {'no', 'gausswin', 'boxcar'});
cfg = ft_checkopt(cfg,'winlen', 'double');
cfg = ft_checkopt(cfg,'gaussvar', 'double');

cfg = ft_checkconfig(cfg, 'allowed', {'channelcmb', 'psth', 'latency', 'colorbar', 'colormap', 'interpolate', 'window', 'winlen', 'gaussvar'});

% determine the corresponding indices of the requested channel combinations
cfg.channelcmb = ft_channelcombination(cfg.channelcmb, jpsth.label);
cmbindx = zeros(size(cfg.channelcmb));
for k=1:size(cfg.channelcmb,1)
  cmbindx(k,1) = strmatch(cfg.channelcmb(k,1), jpsth.label, 'exact');
  cmbindx(k,2) = strmatch(cfg.channelcmb(k,2), jpsth.label, 'exact');
end
nCmbs 	   = size(cmbindx,1);
if nCmbs~=1, error('Currently only supported for a single channel combination'); end

% for convenience create a separate variable
if isfield(jpsth,'jpsth')
  dens = squeeze(jpsth.jpsth(cmbindx(1,1),cmbindx(1,2),:,:)); % density
else
  dens = squeeze(jpsth.shiftpredictor(cmbindx(1,1),cmbindx(1,2),:,:)); % density
end  

% get rid of the NaNs at the borders
isNaN = 1;
k = 1;
while isNaN
  if any(isnan(jpsth.psth(cmbindx,k)))
    disp('deleting NaNs at the borders')
    jpsth.psth(:,k) = [];
    dens(k,:) = [];
    dens(:,k) = [];
    jpsth.time(k) = [];
  else
    isNaN = 0;
  end
end

isNaN = 1;
k = size(jpsth.psth,2);
while isNaN
  if any(isnan(jpsth.psth(cmbindx,k)))
    disp('deleting NaNs at the borders')
    jpsth.psth(:,k) = [];
    dens(k,:) = [];
    dens(:,k) = [];
    jpsth.time(k) = [];
    k = size(jpsth.psth,2);    
  else
    isNaN = 0;
  end
end

% select the time
minTime     = jpsth.time(1);
maxTime     = jpsth.time(end);
if strcmp(cfg.latency,'maxperiod')
  cfg.latency = [minTime maxTime];
elseif strcmp(cfg.latency,'poststim')
  cfg.latency = [0 maxTime];
elseif strcmp(cfg.latency,'prestim')
  cfg.latency = [minTime 0];
end

% check whether the time window fits with the data
if (cfg.latency(1) < minTime), cfg.latency(1) = minTime;
  warning('Correcting begin latency of averaging window');
end
if (cfg.latency(2) > maxTime), cfg.latency(2) = maxTime;
  warning('Correcting end latency of averaging window');
end

% get the samples of our window, and the binwidth of the JPSTH
timeSel       = jpsth.time>=cfg.latency(1) & jpsth.time <= cfg.latency(2);
sampleTime    = mean(diff(jpsth.time)); % get the binwidt


% smooth the jpsth with a kernel if requested
if ~strcmp(cfg.window, 'no')
  
  % construct the kernel
  winTime       = [fliplr(0:-sampleTime:-cfg.winlen) sampleTime:sampleTime:cfg.winlen];
  winLen        = length(winTime);
  if strcmp(cfg.window, 'gausswin') % multivariate gaussian
    A =  winTime'*ones(1,winLen);
    B = A';
    T = [A(:) B(:)]; % makes rows with each time combination on it
    covmat  = diag([cfg.gaussvar cfg.gaussvar]); % covariance matrix
    win    = mvnpdf(T,0,covmat); % multivariate gaussian function
  elseif strcmp(cfg.window, 'boxcar')
    win    = ones(winLen);
  end
  
  % turn into discrete probabilities again (sum(p)=1);
  win = win./sum(win(:));
  win = reshape(win,[],winLen);
  
  % do 2-D convolution and rescale
  dens       = conv2(dens,win,'same');
  outputOnes = conv2(ones(size(dens)),win,'same'); 
  rescale    = 1./outputOnes;
  dens       = dens.*rescale;
end

bins = jpsth.time;
ax(1) = newplot;  % create a new axis object
origPos = get(ax(1), 'Position');
pos     = origPos;
if cfg.interpolate>1
  binAxis = linspace(min(bins), max(bins), round(length(bins)*cfg.interpolate)); % CHANGE THIS
  densInterp = interp2(bins(:), bins(:), dens, binAxis(:), binAxis(:)', 'spline');
  jpsthHdl = imagesc(binAxis, binAxis, densInterp);
else
  jpsthHdl = imagesc(bins,bins,dens);
end
axis xy;

colormap(cfg.colormap); % set the colormap
view(ax(1),2); % use the top view
grid(ax(1),'off');  % toggle grid off

% we need to leave some space for the colorbar on the top
if strcmp(cfg.colorbar,'yes'), pos(3) = pos(3)*0.95; end

% plot the histogram if requested
if strcmp(cfg.psth,'yes')
  
  startPos = pos(1:2) + pos(3:4)*0.2 ; % shift the height and the width 20%
  sz       = pos(3:4)*0.8;             % decrease the size of width and height
  set(ax(1), 'ActivePositionProperty', 'position', 'Position', [startPos sz])
  
  % scale the isi such that it becomes 1/5 of the return plot
  for iChan = 1:2
    % get the data and create the strings that should be the labels
    psth     = jpsth.psth(cmbindx(1,iChan),timeSel);
    label    = jpsth.label{cmbindx(1,iChan)};
    if iChan==2
      startPos = pos(1:2) + pos(3:4).*[0.2 0] ; % shift the width 20%
      sz       = pos(3:4).*[0.8 0.2]; % and decrease the size of width and height
      ax(2) = axes('Units', get(ax(1), 'Units'), 'Position', [startPos sz],...
        'Parent', get(ax(1), 'Parent'));
      psthHdl(iChan) = bar(jpsth.time(timeSel),psth,'k'); % plot the psth under the jpsth
      set(ax(2),'YDir', 'reverse')
      ylabel(label)
    else
      startPos = pos(1:2) + pos(3:4).*[0 0.2] ; % shift the height 20%
      sz       = pos(3:4).*[0.2 0.8]; % and decrease the size of width and height
      ax(3) = axes('Units', get(ax(1), 'Units'), 'Position', [startPos sz],...
        'Parent', get(ax(1), 'Parent'));
      psthHdl(iChan) = barh(jpsth.time(timeSel),psth,'k'); % plot the psth left
      set(ax(3),'XDir', 'reverse')
      xlabel(label)
    end
  end
  set(ax(2), 'YAxisLocation', 'Right', 'XGrid', 'off'); % change the y axis location
  set(ax(3), 'XAxisLocation', 'Top', 'XGrid', 'off');   % change the x axis location
  set(psthHdl, 'BarWidth', 1); % make sure bars have no space between
  
  % change the limits to get the same time limits
  set(ax(2),'XLim', cfg.latency)
  set(ax(3),'YLim', cfg.latency)
  set(get(ax(2),'XLabel'),'String', 'time (sec)')
  set(get(ax(3),'YLabel'),'String', 'time (sec)')
  set(ax(1), 'YTickLabel', {},'XTickLabel', {}); % we remove the labels from JPSTH now
  H.psthleft = psthHdl(2);
  H.psthbottom = psthHdl(1);
  
elseif strcmp(cfg.psth,'no')
  xlabel('time (sec)')
  ylabel('time (sec')
  set(ax(1), 'ActivePositionProperty', 'position', 'Position',pos)
end

% create the colorbar if requested
if strcmp(cfg.colorbar,'yes')
  caxis([min(dens(:))-0.05 max(dens(:))+0.05])
  colormap(cfg.colormap);                  % create the colormap as the user wants
  H.colorbarHdl = colorbar;                % create a colorbar
  
  % create a position vector and reset the position, it should be alligned right of JPSTH
  startPos = [(pos(1) + 0.96*origPos(3)) (pos(2) + pos(4)*0.25)];
  sizePos  = [0.04*origPos(3) 0.75*pos(4)];
  set(H.colorbarHdl, 'Position', [startPos sizePos])
  
  % set the text, try all the possible configurations
  try
    isNormalized =  strcmp(jpsth.cfg.normalization,'yes');
  catch
    isNormalized = 0;
  end
  try
    unit = jpsth.cfg.previous.outputunit;
    if strcmp(unit,'rate')
      isRate = 1;
    elseif strcmp(unit,'spikecount')
      isRate = 2;
    end
  catch
    isRate = 3;
  end
  if isNormalized
    colorbarLabel = 'Normalized Joint Activity';
  elseif isRate==1
    colorbarLabel = 'Joint Firing Rate (spikes^2/sec^2)';
  elseif isRate==2
    colorbarLabel = 'Joint Spike Count (spikes^2)';
  else
    colorbarLabel = 'Joint Peristimulus Activity';
  end
  try
    if strcmp(jpsth.cfg.shiftpredictor,'yes')
      colorbarLabel = strcat(colorbarLabel,' Shift Predictor');
    end
  catch,end
  set(get(H.colorbarHdl,'YLabel'),'String',colorbarLabel)
end

set(ax,'TickDir','out')
set(ax(1),'XLim', cfg.latency,'YLim', cfg.latency)

% get the psth limits for constraining the zooming and panning
psthLim = {};
if strcmp(cfg.psth,'yes')
  psthLim{1}   = get(ax(2),'YLim');
  psthLim{2}   = get(ax(3),'XLim');
end

% collect the handles
H.ax    = ax;
H.jpsth = jpsthHdl;
H.cfg   = cfg;

% constrain the zooming and zoom psth together with the jpsth, remove ticklabels jpsth
set(zoom,'ActionPostCallback',{@mypostcallback,ax,cfg.latency,psthLim});
set(pan,'ActionPostCallback',{@mypostcallback,ax,cfg.latency,psthLim});

% do the general cleanup and bookkeeping at the end of the function

ft_postamble previous jpsth
ft_postamble provenance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = mypostcallback(fig,evd,ax,lim,psthLim)

currentAxes = evd.Axes;
xlim = get(currentAxes, 'XLim');
ylim = get(currentAxes, 'YLim');

if currentAxes==ax(1)
  
  % reset the x limits  to the shared limits
  if lim(1)>xlim(1), xlim(1) = lim(1); end
  if lim(2)<xlim(2), xlim(2) = lim(2); end
  set(ax(1:2), 'XLim',xlim)
  
  % reset the y limits
  if lim(1)>ylim(1), ylim(1) = lim(1); end
  if lim(2)<ylim(2), ylim(2) = lim(2); end
  set(ax([1 3]), 'YLim',ylim)
  
elseif currentAxes == ax(2)
  % reset the x and y limits
  if lim(1)>xlim(1), xlim(1) = lim(1); end
  if lim(2)<xlim(2), xlim(2) = lim(2); end
  
  if psthLim{1}(1)>ylim(1), ylim(1) = psthLim{1}(1); end
  if psthLim{1}(2)<ylim(2), ylim(2) = psthLim{1}(2); end
  
  set(ax(1:2), 'XLim',xlim)
  set(ax(3), 'YLim', xlim)
  set(ax(2),'YLim', ylim)
  
elseif currentAxes == ax(3)
  
  % y limits are shared limits here
  if lim(1)>ylim(1), ylim(1) = lim(1); end
  if lim(2)<ylim(2), ylim(2) = lim(2); end
  
  % psth limit is specific to its x axis
  if psthLim{2}(1)>xlim(1), xlim(1) = psthLim{2}(1); end
  if psthLim{2}(2)<xlim(2), xlim(2) = psthLim{2}(2); end
  
  set(ax(1:2), 'XLim',ylim)
  set(ax(3), 'YLim', ylim)
  set(ax(3),'XLim', xlim)
  
end

set(ax(1), 'YTickLabel', {});
set(ax(1), 'XTickLabel', {});

