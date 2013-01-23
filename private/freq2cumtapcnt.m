function freq = freq2cumtapcnt(freq,fsample)

hasrpt   = ~isempty(strfind(freq.dimord, 'rpt'));
hastim   = ~isempty(strfind(freq.dimord, 'time'));

if ~hasrpt,
 error('computation of number of tapers is not possible, there is not enough information');
end

if hastim,
  tapsmo = freq.cfg.tapsmofrq;
  Nfrq   = length(freq.freq);
  Nsmp   = round(freq.cfg.t_ftimwin.*fsample);
  Ntim   = length(freq.time);
  for j = 1:Nfrq 
    %number of tapers per frequency
    dum       = dpss(Nsmp(j),Nsmp(j).*(tapsmo(j)/fsample));
    numtap(j) = size(dum,2) - 1 ;
  end
  if isfield(freq, 'fourierspctrm') && all(numtap==numtap(1)),
    Ntrl = size(freq.fourierspctrm,1)./numtap(1);
  else
    Ntrl = size(freq.powspctrm,1);
  end
  cumtapcnt = nan(Ntrl, Nfrq, Ntim);    
  for j = 1:Nfrq
    for k = 1:Ntim
      if isfield(freq, 'fourierspctrm'), 
        sel = find(~isnan(freq.fourierspctrm(1:numtap(1):end, 1, j, k)));
      elseif isfield(freq, 'powspctrm'),
        sel = find(~isnan(freq.powspctrm(:, 1, j ,k)));
      end
      cumtapcnt(sel, j, k) = numtap(j);
    end
  end
  freq.cumtapcnt = cumtapcnt;
else
  Nrpt   = size(freq.powspctrm,1);
  tapsmo = freq.cfg.tapsmofrq;
  for j = 1:Nrpt
    dum = dpss(freq.cumsumcnt(j),freq.cumsumcnt(j)*(tapsmo/fsample));
    freq.cumtapcnt(j) = size(dum,2) - 1;
  end
  freq.cumtapcnt = freq.cumtapcnt(:);
end



