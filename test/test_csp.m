function test_csp
% Please beware of notations mistakes in [1]. For example equation (1) does not
% compute the t by t channel covariance, but an n by n time covariance matrix.

% [1] Zoltan J. Koles. The quantitative extraction and topographic mapping of
%     the abnormal components in the clinical EEG. Electroencephalography and
%     Clinical Neurophysiology, 79(6):440--447, December 1991.

clear; close all;
% Create signals with different variance. We use a degenerate covariance
% structure to stress the whitening.
p = 6; n = 100; m = 4;
S1 = diag([0 1 1 1 1 3]) * randn(p, n);
S2 = diag([0 1 1 1 1 .1]) * randn(p, n);

% randomly mix signals
A = randn(p, p);
X1 = A * S1;
X2 = A * S2;

% get covariance
C1 = cov(X1');
C2 = cov(X2');

% find unmixing matrix
W = csp(C1, C2, m);

% test CSP properties
D1 = W * C1 * W';
D2 = W * C2 * W';

assert(all(diff(diag(D1)) <= 0), ...
  'CSP variance is not descending for condition 1.');
assert(norm(D1 + D2 - eye(m)) < 1e-10, ...
  'CSP does not whiten correctly.');
assert(norm(diag(diag(D1)) - D1) < 1e-10, ...
  'CSP does not diagonalize correctly.');
