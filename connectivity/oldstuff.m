%the following has something to do with computing statistics of tfrs of
%pdc/dtf, with respect to a baseline interval

%if ~isfield(cfg, 'cohmethod'), cfg.cohmethod = 'coh';           end;
%if ~iscell(cfg.cohmethod),     cfg.cohmethod = {cfg.cohmethod}; end;
%if ~isfield(cfg, 'submethod'), cfg.submethod = '';              end;
%if ~isempty(cfg.submethod) && ~isfield(cfg, 'baseline'),
%  cfg.baseline = 'all';
%end
%
%if isfield(cfg, 'baseline') && strcmp(cfg.baseline, 'all'),
%  cfg.baseline = [freq.time(1) freq.time(end)];
%end
%
%if isfield(cfg, 'baseline'),
%  baselineindx = [nearest(freq.time, cfg.baseline(1)) nearest(freq.time, cfg.baseline(2))];
%end
%
%
%if hasrpt,
%  nrpt = size(freq.cumtapcnt, 1);
%else
%  nrpt = 1;
%  dum  = zeros([1 size(freq.crsspctrm)]); dum(1,:,:,:,:) = freq.crsspctrm; freq.crsspctrm = dum;
%  dum  = zeros([1 size(freq.powspctrm)]); dum(1,:,:,:,:) = freq.powspctrm; freq.powspctrm = dum;
%  dum  = zeros([1 size(freq.transfer) ]); dum(1,:,:,:,:) = freq.transfer;  freq.transfer  = dum;
%  dum  = zeros([1 size(freq.itransfer)]); dum(1,:,:,:,:) = freq.itransfer; freq.itransfer = dum;
%  dum  = zeros([1 size(freq.noisecov) ]); dum(1,:,:,:,:) = freq.noisecov;  freq.noisecov  = dum;
%  hasrpt = 1;
%end
%if hastim,
%  ntoi = length(freq.time);
%else
%  ntoi = 1;
%end
%nfoi  = length(freq.freq);
%nchan = length(freq.label);
%ncmb  = size(freq.labelcmb,1);
%ntap  = freq.cumtapcnt(1);
%
%for m = 1:length(cfg.cohmethod)
%  switch cfg.cohmethod{m}
%    case {'coh' 'coh2'}
%      for k = 1:ncmb
%        cmbindx(k,1) = match_str(freq.label,freq.labelcmb(k,1));
%        cmbindx(k,2) = match_str(freq.label,freq.labelcmb(k,2));
%      end
%
%      sumcohspctrm = zeros([ncmb  nfoi ntoi]);
%      sumpowspctrm = zeros([nchan nfoi ntoi]);
%      sqrcohspctrm = zeros([ncmb  nfoi ntoi]);
%      sqrpowspctrm = zeros([nchan nfoi ntoi]);
%      warning off;
%      for n = 1:nrpt
%        crsspctrm    = abs(reshape(mean(freq.crsspctrm(n,:,:,:,:),5), [ncmb  nfoi ntoi]));
%        tmppowspctrm = abs(reshape(mean(freq.powspctrm(n,:,:,:,:),5), [nchan nfoi ntoi]));
%
%   if strcmp(cfg.cohmethod{m}, 'coh'),
%     tmpcohspctrm = crsspctrm./sqrt(abs(tmppowspctrm(cmbindx(:,1),:,:,:)).*abs(tmppowspctrm(cmbindx(:,2),:,:,:)));
%        else
%          tmph = reshape(freq.transfer(n,:,:,:,:), [nchan nchan nfoi ntoi ntap]);
%     for flop = 1:nfoi
%       for tlop = 1:ntoi
%         dum                       = tmph(:,:,flop,tlop)*tmph(:,:,flop,tlop)';
%         tmpcohspctrm(:,flop,tlop) = reshape(dum./sqrt(abs(diag(dum))*abs(diag(dum))'), [ncmb 1]);
%       end
%     end
%   end
%
%   if ~isempty(cfg.submethod), tmpcohspctrm = baseline(tmpcohspctrm, cfg.submethod, baselineindx); end
%        if ~isempty(cfg.submethod), tmppowspctrm = baseline(tmppowspctrm, cfg.submethod, baselineindx); end
%   sumcohspctrm = tmpcohspctrm    + sumcohspctrm;
%   sqrcohspctrm = tmpcohspctrm.^2 + sqrcohspctrm;
%   sumpowspctrm = tmppowspctrm    + sumpowspctrm;
%   sqrpowspctrm = tmppowspctrm.^2 + sqrpowspctrm;
%      end
%      warning on;
%      cohspctrm = sumcohspctrm./nrpt;
%      powspctrm = sumpowspctrm./nrpt;
%
%      if nrpt>1,
%        bias         = (nrpt - 1)^2;
%        cohspctrmvar = bias.*(sqrcohspctrm - (sumcohspctrm.^2)/nrpt)./(nrpt-1);
%        powspctrmvar = bias.*(sqrpowspctrm - (sumpowspctrm.^2)/nrpt)./(nrpt-1);
%        cohspctrmsem = sqrt(cohspctrmvar./nrpt);
%        powspctrmsem = sqrt(powspctrmvar./nrpt);
%      end
%    case 'dtf'
%      sumdtf = zeros(ncmb, nfoi, ntoi, ntap);
%      sqrdtf = zeros(ncmb, nfoi, ntoi, ntap);
%      for n = 1:nrpt
%        tmph   = reshape(freq.transfer(n,:,:,:,:), [nchan nchan nfoi ntoi ntap]);
%        den    = sum(abs(tmph).^2,2);
%        tmpdtf = abs(tmph)./sqrt(repmat(den, [1 nchan 1 1 1]));
%        tmpdtf = reshape(tmpdtf, [ncmb nfoi ntoi ntap]);
%        if ~isempty(cfg.submethod), tmpdtf = baseline(tmpdtf, cfg.submethod, baselineindx); end
%        sumdtf = sumdtf + tmpdtf;
%   sqrdtf = sqrdtf + tmpdtf.^2;
%      end
%      dtf = sumdtf./nrpt;
%
%      if nrpt>1,
%        bias   = (nrpt - 1).^2;
%   dtfvar = bias.*(sqrdtf - (sumdtf.^2)/nrpt)./(nrpt-1);
%   dtfsem = sqrt(dtfvar./nrpt);
%      end
%    case 'pdc'
%      sumpdc = zeros(ncmb, nfoi, ntoi, ntap);
%      sqrpdc = zeros(ncmb, nfoi, ntoi, ntap);
%      for n = 1:nrpt
%        invh = reshape(freq.itransfer(n,:,:,:,:), [nchan nchan nfoi ntoi ntap]);
%        %invh = zeros(size(h));
%        %for j = 1:nfoi
%        %  for k = 1:ntoi
%   %    invh(:,:,j,k) = inv(h(:,:,j,k));
%   %  end
%        %end
%        den    = sum(abs(invh).^2,1);
%        tmp    = abs(invh)./sqrt(repmat(den, [nchan 1 1 1 1]));
%        tmppdc = reshape(tmp, [ncmb nfoi ntoi ntap]);
%        if ~isempty(cfg.submethod), tmppdc = baseline(tmppdc, cfg.submethod, baselineindx); end
%   sumpdc = sumpdc + tmppdc;
%   sqrpdc = sqrpdc + tmppdc.^2;
%      end
%      pdc = sumpdc./nrpt;
%
%      if nrpt>1,
%        bias   = (nrpt - 1).^2;
%   pdcvar = bias.*(sqrpdc - (sumpdc.^2)/nrpt)./(nrpt-1);
%   pdcsem = sqrt(pdcvar./nrpt);
%      end
%    otherwise
%      error('unknown cohmethod specified in cfg.cohmethod');
%  end
%end
%
%%---create output-structure
%fd = [];
%fd.label = freq.label;
%fd.labelcmb = freq.labelcmb;
%fd.freq     = freq.freq;
%if hastim, fd.time = freq.time; end
%fd.nobs     = nrpt;
%fd.dimord   = 'chan_freq_time';
%
%try, fd.pdc       = pdc;       end
%try, fd.pdcsem    = pdcsem;    end
%try, fd.dtf       = dtf;       end
%try, fd.dtfsem    = dtfsem;    end
%try, fd.cohspctrm = cohspctrm; end
%try, fd.powspctrm = powspctrm; end
%try, fd.cohspctrmsem = cohspctrmsem; end
%try, fd.powspctrmsem = powspctrmsem; end
%try, cfg.previous    = freq.cfg;     end
%fd.cfg = cfg;
%
%%---subfunction to do baseline correction
%function [output] = baseline(input, method, baseline)
%
%switch method,
%  case 'relchange'
%    b      = mean(input(:,:,baseline(1):baseline(2)),3);
%    output = input./repmat(b, [1 1 size(input,3) 1]) - 1;
%  case 'diff'
%    b      = mean(input(:,:,baseline(1):baseline(2)),3);
%    output = input-repmat(b, [1 1 size(input,3) 1]);
%  otherwise
%    error('specified baseline-method is not yet implemented');
%end
