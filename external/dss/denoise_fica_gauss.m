function [params, s_new] = denoise_fica_gauss(params, s, state)
% FastICA gauss nonlinearity as DSS denoising function
%   [params, s_new] = denoise_fica_gauss(params, s, state)
%     params  Function specific modifiable parameters
%     params.a  Scale constant (default: 1)
%     state   DSS algorithm state
%     s       Source signal estimate, matrix of row vector signals
%     s_new   Denoised signal estimate

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'FastICA gaussian nonlinearity';
    params.description = '';
    params.param = {'a'};
    params.param_value ={1};
    params.param_type = {'scalar'};
    params.param_desc = {'Scaling constant'};
    params.approach = {'defl','symm'};
    params.alpha = {};
    params.beta = {};
    return;
end

if ~isfield(params, 'a')
  params.a = 1;
end
a = params.a;

s2=s.^2;
ex=exp(-a * s2/2);
gauss =  s.*ex;
beta = mean((1 - a * s2) .*ex, 2);
s_new = gauss - diag(beta) * s;
