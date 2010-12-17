function [hdl] = ft_spike_plot_isireturn(cfg,isih)

% FT_SPIKE_PLOT_ISIRETURN makes a return plot from ISIH structure (output from FT_SPIKE_ISIHIST). A
% return plot (or Poincare plots) plots the isi to the next spike versus the isi from the next 
% spike to the second next spike, and thus gives insight in the second order isi statistics.
% This func also plots the raw isi-histogram on left and bottom and thereby give a rather
% complete visualization of the spike-train interval statistics.
%
%   Inputs:
%     ISIH must be the output structure from SPIKE_ISIH and contain the field 
%     ISIH.isi. If cfg.isihist = 'yes', the field ISIH.isih & ISIH.time must be
%     present as well.
%
%   Use must be as follows:
%      [HDL]=SPIKE_ISIRETURNPLOT(CFG,DATA) 
%
%   General configurations (cfg):
%
%   cfg.spikechannel     = string or index of single spike channel to trigger on (default = 1)
%                          Only one spikechannel can be plotted at a time.
%   cfg.density          = 'yes' or 'no', if 'yes', we will use color shading on top of
%                          the individual datapoints to indicate the density.
%   cfg.scatter          = 'yes' (default) or 'no'. If 'yes', we plot the individual values.
%   
%   General configurations related to smoothing the scatterplot
%
%   cfg.smoothmethod     = 'kernel' (default) or 'hist'.
%                           If 'kernel', we overlay a smooth density plot calculated by 
%                           non-parametric kernel smoothing with cfg.kernel.
%                           If 'hist', we overlay a 2-D histogram.
%   cfg.dt               =  resolution of the 2-D histogram, or of the kernel plot. Since we 
%                           have to smooth for a finite number of values, cfg.dt determines
%                           the resolution of our smooth density plot.
%   cfg.colormap         = N-by-3 colormap (see COLORMAP). Default = hot(256);
%   cfg.interpolate      = 'yes' or 'no', determines whether we interpolate the density
%                           plot
%
%   Specific configurations related to kernel smoothing of scatterplot.
%   cfg.kernel           = 'gausswin' or 'boxcar', or N-by-N matrix containing window
%                           values with which we convolve the scatterplot that is binned
%                           with resolution cfg.dt. N should be uneven, so it can be centered
%                           at each point of the lattice.
%                           'gausswin' is N-by-N multivariate gaussian, where the diagonal of the 
%                           covariance matrix is set by cfg.gaussvar.
%                           'boxcar' is N-by-N rectangular window.
%                           If cfg.kernel is numeric, it should be of size N-by-N.
%   cfg.gaussvar         =  variance  (default = 1/16 of window length in sec).
%   cfg.winlen           =  window length in seconds (default = 5*cfg.dt). The total
%                           length of our window is 2*round*(cfg.winlen/cfg.dt) +1;


if nargin~=2, error('MATLAB:spikestation:isireturn:nargin','Two input arguments required'), end

% general configuration defaults
defaults.spikechannel = {1};                   
defaults.scatter      = {'yes' 'no'};               
defaults.density      = {'yes' 'no'};               
defaults.colormap     = flipud(hot(300)); defaults.colormap  = {defaults.colormap(1:256,:)};
defaults.interpolate  = {'yes' 'no'};               
defaults.scattersize  = {0.3};                 
defaults.smoothmethod = {'kernel' 'hist'};            
defaults.dt           = {0.001};               
defaults.kernel       = {'mvgauss'};          
try
  defaults.winlen       = {cfg.dt*5};
catch
  defaults.winlen       = {0.003};
end
try 
  defaults.gaussvar     = {(diff(cfg.winlen)/4).^2}; 
catch
  defaults.gaussvar     = {(defaults.winlen{1}/4).^2}; 
end
cfg = ft_spike_sub_defaultcfg(cfg,defaults);

% check if all the required fields are there
if ~all(isfield(isih,{'isi' 'label' 'time'}))
    error('MATLAB:spikestation:plot_isireturn:cfg:spikechannel:missingFields',...
          'input ISIH should contain the fields isi, label and time')
end

% get the spikechannels: maybe replace this by one function with checking etc. in it
cfg.channel = ft_channelselection(cfg.spikechannel, isih.label);
spikesel    = match_str(isih.label, cfg.channel);
nUnits      = length(cfg.spikechannel); % number of spike channels
if nUnits~=1, error('MATLAB:spikestation:plot_isireturn:cfg:spikechannel:notOneChan',...
                    'Only one unit can be selected at a time'); 
end  
isi = isih.isi{spikesel};
	
% create the axis 
ax(1) = newplot;    
[origPos,pos] = deal(get(ax(1), 'Position'));
if strcmp(cfg.density,'yes'), pos(3) = pos(3)*0.95; end % because of the color bar which takes place

bins       = min(isih.time):cfg.dt:max(isih.time);
nbins      = length(bins);

% two-dimensional kernel smoothing to get a density behind our scatterplot
if strcmp(cfg.density,'yes') 
  isivld1 = (~isnan(isi(1:end-1))&~isnan(isi(2:end)));  
  isivld2 = (isi(1:end-1)<=max(isih.time)-cfg.dt&isi(2:end)<=max(isih.time)-cfg.dt);
  isivld  = isivld1&isivld2;
  
  hold on
  % make a 2-D histogram, we need this for both the kernel and the hist      
  [N,indx1]  = histc(isi(find(isivld)),bins);
  [N,indx2]  = histc(isi(find(isivld)+1),bins);

  % remove the outer counts
  rmv = indx1==length(bins)|indx2==length(bins);
  indx1(rmv) = [];
  indx2(rmv) = [];    
  dens = full(sparse(indx2,indx1,ones(1,length(indx1)),nbins,nbins));
    
  if strcmp(cfg.smoothmethod,'kernel')
    if cfg.winlen<cfg.dt, error('MATLAB:spikestation:plot_isireturn:cfg:dt:winlen',...
      'please configure cfg.winlen such that cfg.winlen>=cfg.dt')
    end
    winTime       = [fliplr(0:-cfg.dt:-cfg.winlen) cfg.dt:cfg.dt:cfg.winlen];
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
        error('MATLAB:spikestation:plot_isireturn:cfg:kernel:wrongInput',...
          'cfg.kernel should be "gausswin", "boxcar" or numerical N-by-N matrix')
    else 
        win    = cfg.kernel;
        szWin  = size(cfg.kernel);
        if  szWin(1)~=szWin(2)||~mod(szWin(1),2), 
            error('MATLAB:spikestation:spike_isireturnplot:cfg:kernel:wrongSize', ...
            'cfg.kernel should be N-by-N matrix with N an uneven number')
        end      
    end

    % turn into discrete probabilities again (sum(p)=1);
    win = win./sum(win);
    win = reshape(win,[],length(winTime)); % reshape to matrix corresponding to grid again

    % do 2-D convolution and rescale
    dens = conv2(dens,win,'same');
    outputOnes = conv2(ones(size(dens)),win,'same');
    rescale = 1./outputOnes;
    dens = dens.*rescale;                
  end

  % create the surface
  hdl.density = imagesc(bins,bins,dens);
  if strcmp(cfg.interpolate,'yes')
    	binAxis = linspace(min(bins)-0.5*(bins(2)-bins(1)), max(bins)+0.5*(bins(2)-bins(1)), 1000);
        densInterp = interp2(bins(:), bins(:), dens, binAxis(:), binAxis(:)', 'spline');
        hdl.density = imagesc(binAxis, binAxis, densInterp);
  end
  if isrealmat(cfg.colormap) && size(cfg.colormap,2)==3
    colormap(cfg.colormap);
  else
    error('MATLAB:spikestation:plot_isireturn:cfg:colormap', ...
    'cfg.colormap should be N-by-3 numerical matrix')
  end
  view(ax(1),2); % use the top view
  grid(ax(1),'off')    
end
hold on
% create scatter return plot and get the handle
if strcmp(cfg.scatter, 'yes') 
  hdl.scatter = plot(isi(1:end-1), isi(2:end),'ko');
  set(hdl.scatter,'MarkerFaceColor', 'k','MarkerSize',  cfg.scattersize)
end

%keyboard
% plot the histogram (user can cut away in adobe of canvas himself if needed)
posisih = [pos(1:2) + pos(3:4)*0.2 pos(3:4)*0.8]; % shift the height and the width 20%
set(ax(1), 'ActivePositionProperty', 'position', 'Position', posisih)        

startPos = pos(1:2) + pos(3:4).*[0.2 0] ; % shift the width 20%
sz       = pos(3:4).*[0.8 0.2]; % and decrease the size of width and height
ax(2) = axes('Units', get(ax(1), 'Units'), 'Position', [startPos sz],...
'Parent', get(ax(1), 'Parent'));

isihist = isih.avg(spikesel,:); %divide by proportion and 5 to get 20% of size
hdl.isi(1) = bar(isih.time(:),isihist,'k');    
set(ax(2),'YDir', 'reverse')

startPos = pos(1:2) + pos(3:4).*[0 0.2] ; % shift the height 20%
sz       = pos(3:4).*[0.2 0.8]; % and decrease the size of width and height
ax(3) = axes('Units', get(ax(1), 'Units'), 'Position', [startPos sz],...
'Parent', get(ax(1), 'Parent'));

hdl.isi(2) = barh(isih.time(:),isihist,'k');    
set(hdl.isi,'BarWidth',1)
set(ax(2:3), 'Box', 'off')
set(ax(3),'XDir', 'reverse')
set(ax(1), 'YTickLabel', {},'XTickLabel',{}, 'XTick',[],'YTick', []);
  
% take care of all the x and y limits at once
set(ax,'Box', 'off','TickDir','out')
set(ax(1),'XLim', [0 max(isih.time)],'YLim',[0 max(isih.time)]);

limIsi = [0 max(isih.avg(:))*1.05];
set(ax(2),'XLim', [0 max(isih.time)],'YLim',limIsi);
set(ax(3),'YLim', [0 max(isih.time)],'XLim',limIsi);
set(ax(2),'YAxisLocation', 'Right')
set(ax(3),'XAxisLocation', 'Top')
set(get(ax(2),'Xlabel'),'String','isi(n) (sec)')
set(get(ax(3),'Ylabel'),'String','isi(n+1) (sec)')

% create the colorbar (who doesn't want one)
caxis([min(dens(:)) max(dens(:))]); 
colormap(cfg.colormap);                  % create the colormap as the user wants  
colorbarHdl = colorbar;                % create a colorbar  

% create a position vector and reset the position, it should be alligned right of isih
startPos = [(pos(1) + 0.96*origPos(3)) (pos(2) + pos(4)*0.25)];
sizePos  = [0.04*origPos(3) 0.75*pos(4)];
set(colorbarHdl, 'Position', [startPos sizePos])  
set(get(colorbarHdl,'YLabel'),'String','Number of spikes per bin');   % set the text
hdl.colorbar = colorbarHdl;
hdl.ax       = ax;
hdl.cfg      = cfg;

% constrain the zooming and zoom psth together with the isih, remove ticklabels isih
set(zoom,'ActionPostCallback',{@mypostcallback,ax,[0 max(isih.time)],limIsi});
set(pan,'ActionPostCallback',{@mypostcallback,ax,[0 max(isih.time)],limIsi});

function [] = mypostcallback(fig,evd,ax,lim,limIsi)

currentAxes = evd.Axes;
indx = find(currentAxes==ax);  

if currentAxes==ax(1)
  % get the x limits and reset them
  xlim = get(ax(1), 'XLim');
  ylim = get(ax(1), 'YLim');
  [axlim] = [min(xlim(1),ylim(1)) max(xlim(2),ylim(2))]; % make the zoom symmetric
  if lim(1)>axlim(1), axlim(1) = lim(1); end
  if lim(2)<axlim(2), axlim(2) = lim(2); end      
  set(ax(1:2), 'XLim',axlim)
  set(ax([1 3]), 'YLim',axlim);  
elseif currentAxes == ax(2)
  xlim = get(ax(2), 'XLim');
  if lim(1)>xlim(1), xlim(1) = lim(1); end
  if lim(2)<xlim(2), xlim(2) = lim(2); end      
  set(ax(1:2), 'XLim',xlim)
  set(ax(3), 'YLim', xlim)
  
  ylim = get(ax(2), 'YLim');
  if limIsi(1)>ylim(1), ylim(1) = limIsi(1); end
  if limIsi(2)<ylim(2), ylim(2) = limIsi(2); end      
  set(ax(2), 'YLim', ylim)
  set(ax(3), 'XLim', ylim)
  
elseif currentAxes == ax(3)
  ylim = get(ax(3), 'YLim');
  if lim(1)>ylim(1), ylim(1) = lim(1); end
  if lim(2)<ylim(2), ylim(2) = lim(2); end      
  set(ax([1 3]), 'YLim',ylim);
  set(ax(2), 'XLim', ylim)
  
  xlim = get(ax(3), 'XLim');
  if limIsi(1)>xlim(1), xlim(1) = limIsi(1); end
  if limIsi(2)<xlim(2), xlim(2) = limIsi(2); end      
  set(ax(3), 'XLim', xlim)
  set(ax(2), 'YLim', xlim)
end











