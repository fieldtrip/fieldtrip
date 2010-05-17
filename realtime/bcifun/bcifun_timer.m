function bcifun_timer(cfg,data)

% BCIFUN_TIMER plots real time versus packet time

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
