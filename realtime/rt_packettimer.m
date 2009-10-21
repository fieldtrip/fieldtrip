function rt_packettimer(cfg)
% RT_PACKETTIMER can be used to time the rate at which data can be processed
%
% Use as
%
%   rt_packettimer(cfg)
%
% with the following configuration options
%   cfg.bcifun    = processing of the data (default = @bcifun_timer) 
%   cfg.npackets  = the number of packets shown in one plot (default=1000)
%                     after reaching the end
%   cfg.saveplot  = if path is specified, first plot is saved (default=[]);
%   cfg.rellim = y limits of subplot 1 (default = [-100 100])
%
% SEE ALSO:
% rt_process
%
% TO DO:
% - jitter in het binnenhalen van de data; scatterplot!
% - triggers sturen en herhalen (loop closen)
% - tijd schatten waarin matlab nog kan processen
%
% Copyright (C) 2009, Marcel van Gerven
%
% $Log: rt_packettimer.m,v $
% Revision 1.8  2009/04/23 12:50:49  marvger
% update of the BCI realtime code
%
% Revision 1.7  2009/02/27 15:00:57  marvger
% update
%
% Revision 1.6  2009/02/10 12:45:27  marvger
% changed help
%
% Revision 1.5  2009/02/10 12:42:04  marvger
% showing absolute delay now. works well with maxblocksize = 1 in
% rt_fileproxy and cfg.blocksize = [0 1] in trialfun_realtime...
%
% Revision 1.4  2009/02/04 14:47:30  marvger
% move from relative to absolute packet timing
%
% Revision 1.3  2009/02/04 14:36:16  marvger
% changed from cputime to tic-toc; save plot added as option
%
% Revision 1.2  2009/02/04 14:04:06  marvger
% changed handling of global variable
%
% Revision 1.1  2009/02/03 20:25:37  marvger
% timing of packets flow in an online setting
%
% Revision 1.2  2009/02/03 10:13:40  marvger
% rt_timer_cpustart is global to force a reset when rt_timer restarts
%
% Revision 1.1  2009/02/03 09:28:37  marvger
% realtime function for timing of packets
%
%

  if ~isfield(cfg,'bcifun'),    cfg.bcifun = @bcifun_timer; end
  if ~isfield(cfg,'npackets'),  cfg.npackets = 10^2; end
  if ~isfield(cfg,'rellim'),    cfg.rellim = [-1 1]; end
  if ~isfield(cfg,'saveplot'),  cfg.saveplot= []; end
    
  close all;
  f1 = figure(); 
  set(f1,'units','normalized','outerposition',[0 0 1 1]);
  
  % reset persistent variables
  cfg.bcifun();
    
  try    
    rt_process(cfg);
  catch    
    fprintf('%s\n',lasterr);
    close;
  end
  
end

function bcifun_timer(cfg,data)
% plots real time versus packet time
%

  persistent startsample;
  persistent olddf;

  persistent idx;

  if nargin == 0
    idx = []; % reset persistent
    return;
  end

  if isempty(idx)
 
    tic;
    
    startsample = data.endsample;
    
    idx       = 1;
    olddf     = 0;
    
    return;
    
  end
  
  if mod(idx,cfg.npackets) == 1
    
    xl = [floor(idx/cfg.npackets)*cfg.npackets+1 (floor(idx/cfg.npackets)+1)*cfg.npackets];

    hold off
    plot(xl,[0 0],'k--');
    xlim(xl);
    ylim(cfg.rellim);
    xlabel('packet number');
    ylabel('delay (sample time - real time)');
    hold on;
    
  end
  
  cursample = data.endsample;
  curtime   = toc;

  smptime = (cursample - startsample)/data.fsample;
  df = smptime - curtime;
  if df < 0
    plot([idx-1 idx],[olddf df],'r');
  else
    plot([idx-1 idx],[olddf df],'g');
  end
  olddf = df;
  title(sprintf('packet timing for %g s blocks in %d channels at sample frequency %d; sample time = %g, real time = %g, delay = %g',...
    data.blocksize/data.fsample,length(data.label),data.fsample,smptime,curtime,smptime - curtime));  

  idx = idx + 1;

  drawnow;

  if idx==cfg.npackets && strcmp(cfg.saveplot,'yes')
    saveas(gcf,cfg.saveplot,'jpg');
  end
  
end
