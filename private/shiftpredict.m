function [prb, cohobs, mcohrnd] = shiftpredict(cfg, dat, datindx, refindx, trltapcnt);

% SHIFTPREDICT implements a shift-predictor for testing significance
% of coherence within a single condition. This function is a subfunction 
% for SOURCESTATISTICS_SHIFTPREDICT and FREQSTATISTICS_SHIFTPREDICT.
%
% cfg.method
% cfg.numrandomization
% cfg.method
% cfg.method
% cfg.loopdim
% cfg.feedback
% cfg.method
% cfg.loopdim
% cfg.correctm
% cfg.tail

% TODO this function should be reimplemented as statfun_shiftpredict for the general statistics framework

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: shiftpredict.m,v $
% Revision 1.8  2009/10/01 12:55:19  jansch
% some changes
%
% Revision 1.7  2006/04/10 12:15:31  roboos
% added cfg options to the documentation
%
% Revision 1.6  2006/01/30 17:13:30  roboos
% made fprintf feedback slightly more instructive about number of tapers
% use single precision for probability for matlab r14 and later
% make shallow copy of input data (instead of deep copy) if possible
%
% Revision 1.5  2005/08/31 07:48:13  jansch
% included the output of the average shift-predicted coherence
%
% Revision 1.4  2005/08/25 14:41:55  jansch
% robert gelooft het wel
%
% Revision 1.3  2005/08/25 10:41:32  roboos
% added support for a third dimension, in which frequencies can be put
% the coherence is computed between reference and target signals, for each frequency separately
%
% Revision 1.2  2005/08/25 08:10:13  jansch
% many small changes, a.o. added support for multiple tapers, added blockshift subfunction
%
% Revision 1.1  2005/08/20 10:12:26  roboos
% new implementations
%

nsgn = size(dat,1);
ntap = size(dat,2); % total number of tapers over all trials
nfrq = size(dat,3);

if nargin<4
  % assume that each trial consists of a single taper only
  ntrl = size(dat,2);
  trltapcnt = ones(1,ntrl);
  fprintf('assuming one taper per trial\n');
else
  % each trial contains multiple tapers
  ntrl = length(trltapcnt);
  trltapcnt = trltapcnt(:)';
  fprintf('number of tapers varies from %d to %d\n', min(trltapcnt), max(trltapcnt));
end

% allocate memory to hold the probabilities
if version('-release')>=14
  prb_pos  = zeros(length(refindx),length(datindx),nfrq, 'single');
  prb_neg  = zeros(length(refindx),length(datindx),nfrq, 'single');
else
  prb_pos  = zeros(length(refindx),length(datindx),nfrq);
  prb_neg  = zeros(length(refindx),length(datindx),nfrq);
end

% compute the power per taper, per trial, and in total
pow_tap = (abs(dat).^2);
pow_trl = zeros(nsgn,ntrl,nfrq);
for i=1:ntrl
  % this is the taper selection if the number of tapers per trial is different
  tapbeg = 1 + sum([0 trltapcnt(1:(i-1))]);
  tapend =     sum([0 trltapcnt(1:(i  ))]);
  % this would be the taper selection if the number of tapers per trial is identical
  % tapbeg = 1 + (i-1)*ntap;
  % tapend =     (i  )*ntap;
  % average the power per trial over the tapers in that trial
  pow_trl(:,i,:) = mean(pow_tap(:,tapbeg:tapend,:),2);
end
pow_tot = reshape(mean(pow_trl,2), nsgn, nfrq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalise the data, also see the COMPUTECOH subfunction
switch cfg.method
case {'abscoh', 'imagcoh', 'absimagcoh', 'atanh', 'atanh_randphase'}
  % normalise, so that the complex conjugate multiplication immediately results in coherence
  if ~all(trltapcnt==trltapcnt(1))
    error('all trials should have the same number of tapers');
  end
  for i=1:nsgn
    for k=1:nfrq
      dat(i,:,k) = dat(i,:,k) ./ sqrt(pow_tot(i,k)*ntrl*trltapcnt(1));
    end
  end

case {'amplcorr', 'absamplcorr'}
  % normalize so that the multiplication immediately results in correlation
  % this uses the amplitude, which is the sqrt of the power per trial
  dat = sqrt(pow_trl);
  % each signal should have zero mean
  fa = reshape(mean(dat,2), nsgn, nfrq);               % mean over trials
  for i=1:ntrl
    dat(:,i,:) = (dat(:,i,:) - fa);
  end
  % each signal should have unit variance
  fs = reshape(sqrt(sum(dat.^2,2)/ntrl), nsgn, nfrq);  % standard deviation over trials
  for i=1:ntrl
    dat(:,i,:) = dat(:,i,:)./fs;
  end
  % also normalize for the number of trials
  dat = dat./sqrt(ntrl);

otherwise
  error('unknown method for shift-predictor')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the random shuffling index vectors
nrnd    = cfg.numrandomization;
randsel = [];
switch cfg.method
case {'abscoh', 'imagcoh', 'absimagcoh', 'atanh'}
  for i=1:nrnd
    % the number of tapers is identical for each trial
    randsel(i,:) = randblockshift(1:sum(trltapcnt), trltapcnt(1));
  end
case {'amplcorr', 'absamplcorr'}
  for i=1:nrnd
    randsel(i,:) = randperm(ntrl);
  end
case {'atanh_randphase'},
  for i=1:nrnd
    rndphs       = exp(j.*2.*pi.*rand(length(unique(cfg.origtrl))));
    randsel(i,:) = rndphs(cfg.origtrl);
  end
end

% select the subset of reference signals
if all(refindx(:)'==1:size(dat,1))
  % only make a shallow copy to save memory
  datref = dat;
else
  % make a deep copy of the selected data
  datref = dat(refindx,:,:);
end

% select the subset of target signals
if all(datindx(:)'==1:size(dat,1))
  % only make a shallow copy to save memory
  dat = dat;
else
  % make a deep copy of the selected data
  dat = dat(datindx,:,:);
end

% compute the observed coherence
cohobs  = computecoh(datref, dat, cfg.method, cfg.loopdim);
mcohrnd = zeros(size(cohobs));
 
progress('init', cfg.feedback, 'Computing shift-predicted coherence');
for i=1:nrnd
  progress(i/nrnd, 'Computing shift-predicted coherence %d/%d\n', i, nrnd);
  % randomize the reference signal and re-compute coherence
  if all(isreal(randsel(i,:))),
    datrnd = datref(:,randsel(i,:),:);
  else
    datrnd = datref.*conj(repmat(randsel(i,:), [size(datref,1) 1 size(datref,3)]));
  end
  cohrnd = computecoh(datrnd, dat, cfg.method, cfg.loopdim);
  mcohrnd = cohrnd + mcohrnd;
  % compare the observed coherence with the randomized one
  if strcmp(cfg.correctm, 'yes')
    prb_pos = prb_pos + (cohobs<max(cohrnd(:)));
    prb_neg = prb_neg + (cohobs>min(cohrnd(:)));
  else
    prb_pos = prb_pos + (cohobs<cohrnd);
    prb_neg = prb_neg + (cohobs>cohrnd);
  end
end
progress('close');
mcohrnd = mcohrnd./nrnd;

if cfg.tail==1
  clear prb_neg  % not needed any more, free some memory
  prb = prb_pos./nrnd;
elseif cfg.tail==-1
  clear prb_pos  % not needed any more, free some memory
  prb = prb_neg./nrnd;
else
  prb_neg = prb_neg./nrnd;
  prb_pos = prb_pos./nrnd;
  % for each observation select the tail that corresponds with the lowest probability
  prb = min(prb_neg, prb_pos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% this assumes that the data is properly normalised
% and assumes that all trials have the same number of tapers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coh = computecoh(datref, dat, method, loopdim)
nref = size(datref,1);
nsgn = size(dat,1);
nfrq = size(dat,3);
coh = zeros(nref,nsgn,nfrq);
switch loopdim
case 3
  % for each frequency, simultaneously compute coherence between all reference and target signals
  for i=1:nfrq
    switch method
    case 'abscoh'
      coh(:,:,i) = abs(     datref(:,:,i) * dat(:,:,i)');
    case 'imagcoh'
      coh(:,:,i) =     imag(datref(:,:,i) * dat(:,:,i)');
    case 'absimagcoh'
      coh(:,:,i) = abs(imag(datref(:,:,i) * dat(:,:,i)'));
    case {'atanh' 'atanh_randphase'}
      coh(:,:,i) = atanh(abs(datref(:,:,i)* dat(:,:,i)'));
    case 'amplcorr'
      coh(:,:,i) =          datref(:,:,i) * dat(:,:,i)';
    case 'absamplcorr'
      coh(:,:,i) = abs(     datref(:,:,i) * dat(:,:,i)');
    otherwise
      error('unsupported method');
    end
  end

case 1
  % for each reference and target signal, simultaneously compute coherence over all frequencies
  for i=1:nref
  for k=1:nsgn
    switch method
    case 'abscoh'
      coh(i,k,:) = abs(     sum(datref(i,:,:) .* conj(dat(k,:,:)), 2));
    case 'imagcoh'
      coh(i,k,:) =     imag(sum(datref(i,:,:) .* conj(dat(k,:,:)), 2));
    case 'absimagcoh'
      coh(i,k,:) = abs(imag(sum(datref(i,:,:) .* conj(dat(k,:,:)), 2)));
    case {'atanh' 'atanh_randphase'}
      coh(i,k,:) = atanh(abs(sum(datref(i,:,:).* conj(dat(k,:,:)), 2)));
    case 'amplcorr'
      coh(i,k,:) =          sum(datref(i,:,:) .* conj(dat(k,:,:)), 2);
    case 'absamplcorr'
      coh(i,k,:) = abs(     sum(datref(i,:,:) .* conj(dat(k,:,:)), 2));
    otherwise
      error('unsupported method');
    end
  end
  end
otherwise
  error('unsupported loopdim');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that shuffles in blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = randblockshift(n, k);
n    = n(:);
nbin = length(n)/k;
n    = reshape(n, [k nbin]);
n    = n(:,randperm(nbin));
out  = n(:);

