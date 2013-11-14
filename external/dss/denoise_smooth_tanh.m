function [params, s_new] = denoise_smooth_tanh(params, s, state)
% Smooth tanh DSS denoising function
%   [params, s_new] = denoise_smooth_tanh(params, s, state)
%     params  Function specific modifiable parameters
%     state   DSS algorithm state
%     s       Source signal estimate, matrix of row vector signals
%     s_new   Denoised signal estimate

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'Smooth tanh';
    params.description = '';
    params.param = {};
    params.param_value ={};
    params.param_type = {};
    params.param_desc = {};
    params.param_type = {};
    params.approach = {'defl','symm'};
    params.alpha = {};
    params.beta = {'beta_global'};
    return;
end

% smoothed variance
[p, var_tot] = denoise_filter(params, s.^2);

if min(var_tot) < 0
  fprintf('[denoise_smooth_tanh.m] Negative smoothed variance, setting to zero.');
  var_tot(var_tot<0) = 0;
end

var_tot_sqrt = sqrt(var_tot);
mask = 1 - tanh(var_tot_sqrt) ./ var_tot_sqrt;
mask = diag(1./(mean(mask, 2))) * mask;
s_new = s .* mask;
