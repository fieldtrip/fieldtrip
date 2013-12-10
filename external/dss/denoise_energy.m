function [params, s_new] = denoise_energy(params, s, state)
% Energy based DSS denoising function
%   [params, s_new] = denoise_energy(params, s, state)
%     params  Function specific modifiable parameters
%     params.usepow     ... (default: 1.3)
%     params.c          ... (default: ?)
%     params.var_noise  Initial noise variance estimate (default: constant)
%     state   DSS algorithm state
%     s       Source signal estimate, matrix of row vector signals
%     s_new   Denoised signal estimate

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'Energy based denoising';
    params.description = '';
    params.param = {'iternoise', 'usepow', 'var_noise'};
    params.param_value ={1, 1.3, 1};
    params.param_type = {'scalar', 'scalar', 'scalar'};
    params.param_desc = {'Number of iterations for noise estimate', '', 'Initial noise variance estimate.'};
    params.approach = {'defl'};
    params.beta = {'beta_global'};
    return;
end

if ~isfield(params, 'initialized')
  % -- Initialize parameters once
  params.initialized = 1;
  if ~isfield(params,'iternoise')
    params.iternoise = 1;
  end
  if ~isfield(params,'usepow')
    params.usepow = 1.3;
  end
  if ~isfield(params,'gamma')
    [p, var_smooth] = denoise_filter(params, randn(1, length(s)).^2, state);
    noise_var = est_noise_var(var_smooth, 0);
    params.c = 1 / noise_var;
    %fprintf('Noise variance: %f', noise_var);
  end
  if ~isfield(params,'var_noise')
    % Initial value for noise variance estimate is 1
    params.var_noise = 1;
  end
end

[p, var_totsm] = denoise_filter(params, s.^2, state);

for i = 1 : params.iternoise
  params.var_noise = est_noise_var(var_totsm, params.var_noise)*params.c;
end
noise_stretched = repmat(params.var_noise, 1, size(var_totsm, 2));
var_sig = sqrt(noise_stretched.^2 + var_totsm.^2) - noise_stretched;
mask = var_sig.^params.usepow ./ var_totsm;

% TODO: is saturation necessary? (test if ever above 4)
% saturation -> 4
%mask = mask .* (mask < 4) + 4 * (mask >= 4);
if max(mask)>5
  fprintf('[denoise_energy.m] Denoise energy mask saturation: %d\n', max(mask));
end

%mask = mask ./ repmat(mean(mask, 2), 1, size(mask, 2));

s_new = mask .* s;

% --------
function noise_var = est_noise_var(var_tot, var_noise)

noise_var = exp(mean(log(var_tot+repmat(var_noise,1,size(var_tot,2))),2))-var_noise;
