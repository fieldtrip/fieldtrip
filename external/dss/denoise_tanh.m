function [params, s_new] = denoise_tanh(params, s, state)
% Tanh DSS denoising function
%   [params, s_new] = denoise_tanh(params, s, state)
%     params  Function specific modifiable parameters
%     state   DSS algorithm state
%     s       Source signal estimate, matrix of row vector signals
%     s_new   Denoised signal estimate
%
%   Calculates hyperbolic tangent of the signal. Use together with
%   suitable tanh beta function.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'Supergaussian';
    params.description = '';
    params.param = {};
    params.param_value ={};
    params.param_desc = {};
    params.param_type = {};
    params.alpha = {};
    params.beta = {'beta_global', 'beta_tanh'};
    params.approach = {'defl', 'symm'};
    return;
end

s_new = s - tanh(s);
