function [spec, se, wt, varargout] = adaptspec_dpss(yk, lambda, adaptflag, varargin)

% ADAPTSPEC_DPSS Adaptive multitaper spectrum estimate
%
%   [spec, se, wt] = adaptspec_dpss(yk, sk, lambda, adaptflag)
%
% Inputs
%   yk       : [nchan x nfft x ntap] complex, tapered data Fourier transforms
%   lambda   : [ntap x 1] real, eigenvalues of DPSS tapers
%   adaptflag : integer, mode (default=0)
%            0 - unweighted average (all tapers equal)
%            1 - weighted by eigenvalues
%            2 - adaptive multitaper (Thomson)
%            3 - quadratic inverse based curvature correction (Prieto),
%            also provides an estimate of the local slope and quadratic
%            curvature
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
% - The quadratic inverse estimation routine deviates a little bit from the
% original code, in that it uses an adjusted fitting procedure for
% frequency bins that are close to DC and Nyquist

if nargin < 3
  adaptflag = 0;
end
yk                  = permute(yk, [2 1 3]);
[nfft, nchan, ntap] = size(yk);
sk                  = abs(yk).^2;

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
  % adaptive scheme as per Thomson, is in a subfunction, because it is also required for adaptflag==3
  [spec, se, wt] = adaptive(sk, lambda, nfft, nchan, ntap);
elseif adaptflag==3
  % quadratic inverse based local curvature correction of the spectrum, requires the Slepian tapers 
  % time courses, and the weights to be computed as per the adaptflag==2 option. It also returns an 
  % estimate of the slope and local curvature.
  tap = varargin{1};
  nw  = varargin{2};
  [spec, se, wt]        = adaptive(sk, lambda, nfft, nchan, ntap);
  [qispec, slope, quad] = qiinv(yk, wt, tap.', lambda, nw);
  qispec = qispec .';
  slope  = slope .';
  quad   = quad .';
end

% permute back
wt   = ipermute(wt, [2 1 3]);
spec = spec.';
se   = se.';
if adaptflag==3
  varargout{1} = qispec;
  varargout{2} = slope;
  varargout{3} = quad;
end

function [spec, se, wt] = adaptive(sk, lambda, nfft, nchan, ntap)

% frequency sampling axis (unit sample rate assumed)
df = 1 / (nfft - 1);

% variance of sk and avg variance
varsk  = sum(sk, 1) * df;     % [1 x nchan x ntap]
dvar   = mean(varsk, 3);

lambda = shiftdim(lambda(:).', -1); % ensure row vector
bk     = repmat(dvar, [1 1 ntap]) .* (1 - repmat(lambda, [1 nchan 1])); % Thomson Eq 5.1b

% initialize
sbar = (sk(:,:,1) + sk(:,:,2)) / 2; % initial guess
spec = sbar;

rerr  = 1e-09;
mloop = 1000;

for i = 1:mloop
  slast = sbar;

  wt1 = sbar.*sqrt(lambda); %should yield nfft x nchan x ntap matrix, not sure whether robust for old matlab
  wt2 = sbar.*lambda + bk;

  wt = min(wt1 ./ wt2, 1.0);
  
  wtsq   = wt.^2;
  skw    = wtsq .* sk;
  wtsum  = sum(wtsq, 3);
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

function se = wt2dof(wt)
se = 2 * (sum(wt, 3).^2) ./ sum(wt.^2, 3);
