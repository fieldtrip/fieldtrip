function [params, s_new] = denoise_fica_tanh(params, s, state)
% FastICA tanh nonlinearity as DSS denoising function
%   [params, s_new] = denoise_fica_tanh(params, s, state)
%     params  Function specific modifiable parameters
%       a     Scaling constant (default: 1)
%     state   DSS algorithm state
%     s       Source signal estimate, matrix of row vector signals
%     s_new   Denoised signal estimate

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'FastICA tanh nonlinearity';
    params.description = '';
    params.param = {'a'};
    params.param_value ={1};
    params.param_type = {'scalar'};
    params.param_desc = {'Steepness modifier'};
    params.approach = {'defl','symm'};
    params.alpha = {};
    params.beta = {};
    return;
end

if ~isfield(params, 'a')
  params.a = 1;
end
a = params.a;

hypTan = tanh(a * s);
s_new = hypTan - a * repmat(mean(1 - hypTan .^ 2, 2), 1, length(s)) .* s;
