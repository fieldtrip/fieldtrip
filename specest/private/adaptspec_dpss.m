function [spec, se, wt] = adaptspec_dpss(yk, lambda, adaptflag)

% ADAPTSPEC_DPSS Adaptive multitaper spectrum estimate
%
%   [spec, se, wt] = adaptspec_dpss(yk, sk, lambda, adaptflag)
%
% Inputs
%   yk     : [nchan x nfft x ntap] complex, tapered data Fourier transforms
%   lambda   : [ntap x 1] real, eigenvalues of DPSS tapers
%   adaptflag : integer, mode (default=0)
%            0 - unweighted average (all tapers equal)
%            1 - weighted by eigenvalues
%            2 - adaptive multitaper (Thomson)
%
% Outputs
%   spec : [nchan x nfft] spectrum estimate
%   se   : [1 x nfft] effective degrees of freedom
%   wt   : [nchan x nfft x ntap] weights applied to tapers
%
% Notes
% - Follows German Prieto’s code (https://github.com/gaprieto/multitaper).
% - The adaptive scheme iterates until convergence or max mloop.
% - ChatGPT provided a rough translation from Python to MATLAB, code
% optimized for MATLAB + multichannel data by JMS

if nargin < 3
  adaptflag = 0;
end

yk = permute(yk, [2 1 3]);

[nfft, nchan, ntap] = size(yk);
sk                   = abs(yk).^2;
if adaptflag==0
  % average across tapers
  wt     = ones(nfft, nchan, ntap);
  sbar   = mean(sk, 3);
  spec   = sbar;
  se     = 2 * ntap * ones(nfft,1);
elseif adaptflag==1
  % weigh the tapered estimates by the tapers' concentration eigenvalues
  wt     = repmat(shiftdim(lambda(:).', -1), nfft, nchan, 1);
  skwsum = sum(sk.*(wt.^2), 3);
  sbar   = skwsum ./ sum(wt.^2, 3);
  spec   = sbar;
  se     = wt2dof(wt);
elseif adaptflag==2
  
  % Freq sampling (unit sample rate assumed)
  df = 1 / (nfft - 1);

  % Variance of Sk's and avg variance
  varsk  = sum(sk, 1) * df;     % [1 x ntap]
  dvar   = mean(varsk, 3);
  
  lambda = shiftdim(lambda(:).', -1); % ensure row vector
  bk     = repmat(dvar, [1 1 ntap]) .* (1 - repmat(lambda, [1 nchan 1])); % Thomson Eq 5.1b
  
  % Initialize
  sbar = (sk(:,:,1) + sk(:,:,2)) / 2; % initial guess
  spec = sbar;

  rerr  = 1e-12;
  mloop = 1000;
   
  onevec = ones(nfft, 1);
  for i = 1:mloop
    slast = sbar;
    
    for m = 1:ntap
      wt1(:,:,m) = sbar*sqrt(lambda(m)); % [nfft x 1] x [1 x ntap]
      wt2(:,:,m) = sbar*lambda(m) + bk(onevec,:,m);
    end
    wt = min(wt1 ./ wt2, 1.0);

    skw    = (wt.^2) .* sk;
    wtsum  = sum(wt.^2, 3);
    skwsum = sum(skw, 3);
    sbar   = skwsum ./ wtsum;
    oerr   = max(max(abs((sbar - slast) ./ (sbar + slast))));
  
    if i == mloop
      spec = sbar;
      warning('adaptspec did not converge, rerr=%g (target %g)', oerr, rerr);
      break;
    end

    if oerr <= rerr
      spec = sbar;
      break;
    end
  end
  se = wt2dof(wt);
end

% permute back
wt   = ipermute(wt, [2 1 3]);
spec = spec.';
se   = se.';

function se = wt2dof(wt)

se = 2 * (sum(wt, 3).^2) ./ sum(wt.^2, 3);
