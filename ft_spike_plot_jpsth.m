function [H] = plot_jpsth(cfg,jpsth)

% plot_jpsth makes a return plot from JPSTH structure (output from
% plot_jpsth). 
%
%   Inputs:
%     JPSTH must be the output structure from SPIKE_JPSTH and contain the field 
%     JPSTH.avg. If cfg.psth = 'yes', the field JPSTH.psth must be
%     present as well.
%
%   Use must be as follows:
%      [HDL]=plot_jpsth(CFG,JPSTH) 
%
%   Here we should definitely think of making a matrix
%   With the same electrode somehow as the diagonal
%
%   General configurations (cfg):
%
%   cfg.channelcmb       = string or index of single channel combination to trigger on.
%                          See SPIKESTATION_FT_SUB_CHANNELCOMBINATION for details.
%   cfg.psth             = 'yes' (default) or 'no'. If 'yes', the psth histogram is down
%                          the x and left from the y axis.
%   cfg.axlim            = [begin end] in seconds or 'max' (default), 'prestim' or
%                          'poststim';
%   cfg.colorbar         = 'yes' (default) or 'no'   
%   cfg.colormap         =  N-by-3 colormap (see COLORMAP). or 'auto' (default,hot(256))
%   cfg.interpolate      = 'yes' or 'no', determines whether we interpolate the density
%                           plot
%   Configurations related to smoothing 
%
%   cfg.smooth           = 'yes' or 'no' (default) 
%                           If 'yes', we overlay a smooth density plot calculated by 
%                           non-parametric symmetric kernel smoothing with cfg.kernel.
%   cfg.kernel           = 'gausswin' or 'boxcar', or N-by-N matrix containing window
%                           values with which we convolve the scatterplot that is binned
%                           with resolution cfg.dt. N should be uneven, so it can be centered
%                           at each point of the lattice.
%                           'gausswin' is N-by-N multivariate gaussian, where the diagonal of the 
%                           covariance matrix is set by cfg.gaussvar.
%                           'boxcar' is N-by-N rectangular window.
%                           If cfg.kernel is numeric, it should be of size N-by-N.
%   cfg.gaussvar         =  variance  (default = 1/16 of window length in sec).
%   cfg.winlen           =  window length in seconds (default = 5*binwidth). The total
%                           length of our window is 2*round*(cfg.winlen/binwidth)
%                           where binwidth is the binwidth of the jpsth (jpsth.time(2) -
%                           jpsth.time(1)).
%
%   See also SPIKE_JOINTPSTH, SPIKE_PSTH, SPIKE_XCORR
%

if nargin~=2, error('MATLAB:spikestation:plot_jpsth:narginWrong',...
    'Two input arguments are required')
end

% general configuration options
defaults.channelcmb   = {{1,2}};
defaults.psth         = {'yes' 'no'};                            
defaults.latency      = {'maxperiod'};                         
defaults.colorbar     = {'yes' 'no'};                            
defaults.colormap     = {jet(256)};                         
defaults.interpolate  = {'no' 'yes'};                            

% further kernel smoothing configurations
defaults.smooth       = {'no' 'yes'};
defaults.kernel       = {'mvgauss'};
try 
  defaults.winlen       = {5*(jpsth.time(2)-jpsth.time(1))};
catch
  defaults.winlen      = {0.01}; %sec
end
try 
  defaults.gaussvar     = {(cfg.winlen/4).^2}; 
catch
  defaults.gaussvar     = {(defaults.winlen{1}/4).^2}; 
end
cfg = ft_spike_sub_defaultcfg(cfg,defaults);

% channel combination selection
cmbindx = ft_spikestation_sub_channelcombination(cfg.channelcmb, jpsth.label);
nCmbs = size(cmbindx,1);
if nCmbs~=1, error('MATLAB:spikestation:plot_jpsth:cfg:channelcmb:tooManySelected', ...
    'Currently only supported for a single channel combination')
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
elseif ~isrealvec(cfg.latency)||length(cfg.latency)~=2
  error('MATLAB:spikestation:plot_jpsth:cfg:latency',...
    'cfg.latency should be "max" or 1-by-2 numerical vector');
end
if cfg.latency(1)>=cfg.latency(2),
     error('MATLAB:spikestation:plot_jpsth:cfg:latency:wrongOrder',...
    'cfg.latency(2) should be greater than cfg.latency(1)')
end
% check whether the time window fits with the data
if (cfg.latency(1) < minTime), cfg.latency(1) = minTime; 
  warning('MATLAB:spikestation:plot_jpsth:correctLatencyBeg',...
          'Correcting begin latency of averaging window');
end
if (cfg.latency(2) > maxTime), cfg.latency(2) = maxTime;
  warning('MATLAB:spikestation:plot_jpsth:correctLatencyEnd',...
          'Correcting end latency of averaging window');
end

% get the samples of our window, and the binwidth of the JPSTH
timeSel = jpsth.time>=cfg.latency(1) & jpsth.time <= cfg.latency(2);
sampleTime    = jpsth.time(2) - jpsth.time(1); % get the binwidth

% for convenience create a separate variable
dens = squeeze(jpsth.avg(:,:,cmbindx(1,1),cmbindx(1,2))); % density

% smooth the jpsth with a kernel if requested
if strcmp(cfg.smooth,'yes')

  % construct the kernel
  winTime       = [fliplr(0:-sampleTime:-cfg.winlen) sampleTime:sampleTime:cfg.winlen];
  winLen        = length(winTime);    
  if strcmp(cfg.kernel, 'mvgauss') % multivariate gaussian
      A =  winTime'*ones(1,winLen);
      B = A';
      T = [A(:) B(:)]; % makes rows with each time combination on it
      covmat  = diag([cfg.gaussvar cfg.gaussvar]); % covariance matrix
      win    = mvnpdf(T,0,covmat); % multivariate gaussian function
  elseif strcmp(cfg.kernel, 'boxcar')
      win    = ones(winLen);
  elseif ~isrealmat(cfg.kernel)
      error('MATLAB:spikestation:plot_jpsth:cfg:kernel:unknownOption',...
        'cfg.kernel should be "gausswin", "boxcar" or numerical N-by-N matrix')
  else 
      win    = cfg.kernel;
      szWin  = size(cfg.kernel);
      if  szWin(1)~=szWin(2)||~mod(szWin(1),2), 
          error('MATLAB:spikestation:plot_jpsth:cfg:kernel:wrongSize', ...
          'cfg.kernel should be 2-dimensional, N-by-N matrix with N uneven')
      end      
  end
 
  % turn into discrete probabilities again (sum(p)=1);
  win = win./sum(win(:));
  win = reshape(win,[],winLen); 

  % do 2-D convolution and rescale
  dens  = conv2(dens,win,'same');
  outputOnes = conv2(ones(size(dens)),win,'same'); % NOTE $$$ SHOULD BE IN RIGHT ROW/COLUMN!
  rescale = 1./outputOnes;
  dens  = dens.*rescale;                
end

bins = jpsth.time;
ax(1) = newplot;  % create a new axis object  
origPos = get(ax(1), 'Position');
pos     = origPos;
if strcmp(cfg.interpolate,'yes')
    binAxis = linspace(min(bins), max(bins), 1000); % CHANGE THIS
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
  caxis([min(dens(:)) max(dens(:))])  
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

