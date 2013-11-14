function [params, s_new] = denoise_pow3(params, s, state)
% Kurtosis based denoising function
%   [params, s_new] = denoise_pow3(params, s, state)
%     params  Function specific modifiable parameters
%     state   DSS algorithm state
%     s       Source signal estimate, matrix of row vector signals
%     s_new   Denoised signal estimate
%
%   Calculates third power of the signal. Equals kurtosis if source
%   signal has unit variance.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'Kurtosis based denoising';
    params.description = '';
    params.param = {};
    params.param_value ={};
    params.param_type = {};
    params.param_desc = {};
    params.param_type = {};
    params.approach = {'defl','symm'};
    params.alpha = {};
    params.beta = {'beta_global','beta_pow3'};
    return;
end

s_new = s.^3;
