function [mask] = estimate_mask(s, filter_params, iterations)
% Creates binary mask based on SNR estimate of the signal
%   [mask] = estimate_mask(s, filter_params, iterations)
%     s              Source signal
%     filter_params  Filter parameters for smoothing variance estimate
%     iterations     Number of iterations for estimating noise variance

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

filter_h = @denoise_filter;

if nargin<2;
    % default filtering is lowpass DCT
    T = length(s);
    t = 1:T;
    filter_params.filter_dct = exp(-0.5*(t-1).^2 / (T/64).^2 );
end
if nargin<3; iterations = 8; end

% noise variance estimate for gaussian noise
[p, var_smooth] = feval(filter_h, filter_params, randn(1, length(s)).^2, []);
noise_var = est_noise_var(var_smooth, 0);
normalization_c = 1 / noise_var;

% smoothed variance
[p, var_totsm] = feval(filter_h, filter_params, s.^2, []);

% iterate signal noise variance estimate
var_noise = 1;
for i = 1 : iterations
  var_noise = est_noise_var(var_totsm, var_noise)*normalization_c;
end

% create binary mask
mask = var_totsm>var_noise;

%DEBUG
%fprintf('Noise: %d\n', var_noise);
%clf
%subplot(3, 1, 1);
%plot(s);
%subplot(3, 1, 2);
%plot(var_totsm);
%subplot(3, 1, 3);
%plot(mask);
%axis([0 length(s) -0.5 1.5])

% --------
function noise_var = est_noise_var(var_tot, var_noise)
% Estimates noise variance based on total variance estimate and
% previous noise variance estimate.
noise_var = exp(mean(log(var_tot+repmat(var_noise,1,size(var_tot,2))),2))-var_noise;
