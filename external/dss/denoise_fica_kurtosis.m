function [params, s_new] = denoise_fica_kurtosis(params, s, state)
% FastICA kurtosis nonlinearity as DSS denoising function
%   [params, s_new] = denoise_fica_kurtosis(params, s, state)
%     params  Function specific modifiable parameters
%     state   DSS algorithm state
%     s       Source signal estimate, matrix of row vector signals
%     s_new   Denoised signal estimate

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'FastICA kurtosis nonlinearity';
    params.description = '';
    params.param = {};
    params.param_value ={};
    params.param_type = {};
    params.param_desc = {};
    params.approach = {'defl','symm'};
    params.alpha = {};
    params.beta = {};
    return;
end

s_new = s.^3 - 3 * s;
