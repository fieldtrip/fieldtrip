function [qispec, slope, quad] = qiinv(yk, wt, vn, lambda, nw)
% QIINV  Quadratic-spectrum (QI) estimate following Prieto et al. (2007),
% adjusted for the frequency bins within NW from DC and Nyquist
%
%   [qispec, slope, quad] = qiinv(yk, wt, vn, lambda, nw)
%
% Inputs
%   yk    : [nfft x nchan x ntap] complex, single taper Fourier spectra
%   wt    : [nfft x nchan x ntap] real, adaptive weights
%   vn    : [nsmp x ntap] real, Slepian sequences
%   lambda  : [ntap x 1] real, Slepian eigenvalues
%   nw    : scalar, time-bandwidth product
%
% Outputs
%   qispec : [nfft x nchan] QI spectrum estimate (bias-reduced)
%   slope  : [nfft x nchan] first derivative estimate
%   quad   : [nfft x nchan] second derivative estimate
%
% Notes
% - Uses Chebyshev polynomials (unitless basis) as in the reference.
% - Follows the 2021 modification: invert constant term with NNLS first,
%   then invert for derivatives to keep the 2nd derivative independent.
%
% Reference:
%   Prieto, Parker, Thomson, Vernon, Graham (2007), GJI 171, 1269–1281.
%   doi:10.1111/j.1365-246X.2007.03592.x
%
% Authors: MATLAB translation of Python code, by ChatGPT with suggestions from
% Mistral, to adjust the estimates close to DC and Nyquist, checked by JMS

[nsmp, ntap]           = size(vn);
[nfft, nchan, ntap_yk] = size(yk);
assert(ntap_yk==ntap, 'yk and vn must have the same number of tapers (columns).');

nxi   = 79;
L     = ntap^2;

if min(lambda) < 0.9
  warning('careful, poor leakage of eigenvalue %g. ntap may be too large.', min(lambda));
end


% inner bandwidth grid
Vj  = complex(zeros(nxi,  ntap));
bp  = nw / nsmp;                 % W bandwidth
xi  = linspace(-bp, bp, nxi).';  % column vector [nxi x 1]
dxi = xi(2) - xi(1);

% complex valued frequency transform of tapers within the bandwith
for k = 1:ntap
  for i = 1:nxi
    om = 2.0*pi*xi(i);
    [ct, st] = sft(vn(:,k), om); % scalar cosine/sine transforms
    Vj(i,k) = (1.0/sqrt(lambda(k))) * complex(ct, st);
  end
end

% ---------- Build vectorized C and Pk ----------
C  = complex(zeros(L, nchan, nfft));
Pk = complex(zeros(L, nxi));

m = 0;
xk = yk.*wt;
for i = 1:ntap
  for k = 1:ntap
    m = m + 1;
    C(m, :, :)  = transpose(conj(xk(:,:,i)) .* xk(:,:,k)); % [1 x nfft] -> C = ntap^2 * nchan * nfft
    Pk(m, :)    = conj(Vj(:,i)) .* Vj(:,k);                % [1 x nxi]  -> Pk = ntap^2 * gridded bandwidth
  end
end

% trapezoid end weights
Pk(:,1)   = 0.5 * Pk(:,1);
Pk(:,end) = 0.5 * Pk(:,end);

% Chebyshev basis (unitless)
hconstant   = ones(nxi, 1);
hslope = xi / bp;
hquad  = 2.0*( (xi/bp).^2 ) - 1.0;

h1 = (Pk * hconstant) * dxi; % [L x 1]
h2 = (Pk * hslope)    * dxi; % [L x 1]
h3 = (Pk * hquad)     * dxi; % [L x 1]

% hk0 works fine when the tested centre frequency is sufficiently far away from DC or Nyquist
hk0 = [h1, h2, h3]; % [L x 3]
nh  = size(hk0,2);

% QR factorization & covariance of coefficients
[Q0, R0] = qr(hk0);      % NO economy QR: hk = Q*R
Qt0      = Q0';
Leye     = eye(L);
Ri0      = R0 \ Leye;        % same as inv(R), but stable
covb0    = real(Ri0 * Ri0.'); % covariance up to sigma^2 scaling

% preallocate variables
[cte, cte2, slope, quad]               = deal(zeros(nfft, nchan));
[sigma2, cte_var, slope_var, quad_var] = deal(zeros(nfft, nchan));

faxis = (0:(nfft-1))./nfft;
faxis(faxis>0.5) = -1 + faxis(faxis>0.5);

% loop over frequencies and perform regression
cte_out = zeros(1, nchan);
for ii = 1:nfft
  cjk = C(:, :, ii); % [L x nchan]

  f0 = abs(faxis(ii));
  if f0 <= bp || f0 >= 0.5-bp
    [hk, Qt, R, covb] = folded_basis(Pk, xi, f0, bp);
  else
    Qt = Qt0;
    R  = R0;
    covb = covb0;
  end

  % --- Invert constant term with non-negative least squares (NNLS) ---
  % Equivalent to optim.nnls(np.real(h1), np.real(cjk))
  for i = 1:nchan
    cte_out(i) = lsqnonneg(real(h1), real(cjk(:,i)));
  end
  cte2(ii,:) = real(cte_out);
  pred_cte   = h1 * cte2(ii,:);
  cjk2       = cjk - pred_cte;

  % solve for the derivatives given the residual
  btilde = Qt * cjk2;
  hmodel = R \ btilde;         % least-squares solution
  cte(ii,:)   = real(hmodel(1,:));
  slope(ii,:) = real(hmodel(2,:));
  quad(ii,:)  = real(hmodel(3,:));

  pred   = hk * real(hmodel);
  sigma2(ii,:) = sum(abs(cjk2 - pred).^2) / (L - nh);

  cte_var(ii,:)   = sigma2(ii,:) * covb(1,1);
  slope_var(ii,:) = sigma2(ii,:) * covb(2,2);
  quad_var(ii,:)  = sigma2(ii,:) * covb(3,3);
end

% normalize derivatives to physical units
slope     = slope / bp;
quad      = quad  / (bp^2);
slope_var = slope_var / (bp^2);
quad_var  = quad_var  / (bp^4);

% compute the bias and correct
[qispec, qicorr] = deal(zeros(nfft,nchan));

% shrinkage factors, one per term
sr2 = quad.^2 ./ (quad.^2 + quad_var);     % curvature shrinkage
sr1 = slope.^2 ./ (slope.^2 + slope_var);  % slope shrinkage
for ii = 1:nfft
  f0 = abs(faxis(ii));

  if (f0 >= bp) && (f0 <= 0.5 - bp)
    % original correction, as per the paper
    qicorr(ii,:) = sr2(ii,:) .* quad(ii,:) * (1/6) * bp^2;
  else
    % adjusted correction, as per Mistral, correcting for being close to the boundary of the spectrum
    [mu1, mu2] = folded_moments(Pk, xi, f0);
    qicorr(ii,:) = sr1(ii,:) .* slope(ii,:) * mu1 + sr2(ii,:) .* quad(ii,:) * 0.5 * mu2;
  end
  qispec(ii,:) = cte2(ii,:) - qicorr(ii,:);
end

% ===========================
% Helper: single-frequency transform of a real vector v at angular freq om
% Returns real and imag parts separately as cosine and sine sums.
% Matches the Python 'sft' signature: [ct, st] = sft(v, om)
function [ct, st] = sft(v, om)
n = numel(v);
t = (0:n-1).';
c = cos(om .* t);
s = sin(om .* t);
ct = sum(v .* c);      % real (cosine) component
st = sum(v .* s);      % imaginary (sine) component

function [hk, Qt, R, covb] = folded_basis(Pk, xi, f0, bp)

% subfunction to create a slightly adjusted set of basis functions for the
% slope and curvature fits, which is needed for frequency bins that are
% close to DC and Nyquist
dxi = xi(2) - xi(1);

nu = abs(f0 + xi);
nu = 0.5 - abs(0.5 - nu);

hconstant = ones(numel(xi), 1);
hslope = (nu - f0) / bp;
hquad  = 2.0*( (hslope).^2 ) - 1.0;

hk(:,1) = (Pk * hconstant) * dxi; % [L x 1]
hk(:,2) = (Pk * hslope)    * dxi; % [L x 1]
hk(:,3) = (Pk * hquad)     * dxi;

[Q, R]  = qr(hk);
Qt      = Q';
Leye    = eye(size(Q,1));
Ri      = R \ Leye;        % same as inv(R), but stable
covb    = real(Ri * Ri.'); % covariance up to sigma^2 scaling

function [mu1, mu2] = folded_moments(Pk, xi, f0)

% Kernel-weighted moments of (nu - f0) for the folded model.
% Pk: [L x nxi] kernel matrix (complex), xi: [1 x nxi] offsets [-bp, bp]

w  = mean(real(Pk), 1)';  % row-averaged kernel weight [1 x nxi]
nu = abs(f0 + xi);        % DC fold
nu = 0.5 - abs(0.5 - nu); % Nyquist fold
d  = nu - f0;             % model argument minus f0, in f units
W0 = trapz(xi, w);
mu1 = trapz(xi, w .* d)     / W0;
mu2 = trapz(xi, w .* d.^2)  / W0;
