function ft_realtime_packettimer(cfg)

% FT_REALTIME_PACKETTIMER can be used to time the rate at which data can be processed
%
% Use as
%   ft_realtime_packettimer(cfg)
% with the following configuration options
%   cfg.bcifun    = processing of the data (default = @bcifun_timer)
%   cfg.npackets  = the number of packets shown in one plot (default=1000)
%                     after reaching the end
%   cfg.saveplot  = if path is specified, first plot is saved (default=[]);
%   cfg.rellim = y limits of subplot 1 (default = [-100 100])
%
% SEE ALSO:
%   FT_REALTIME_PROCESS

% TO DO:
%   jitter in het binnenhalen van de data; scatterplot!
%   triggers sturen en herhalen (loop closen)
%   tijd schatten waarin matlab nog kan processen

% Copyright (C) 2009, Marcel van Gerven
%
% $Id$

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
  ft_realtime_process(cfg);
catch
  fprintf('%s\n',lasterr);
  close;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCIFUN_TIMER plots real time versus packet time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bcifun_timer(cfg,data)

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

