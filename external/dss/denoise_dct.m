function [params, s_new] = denoise_dct(params, s, state)
% DCT filter function.
%   [s_new, params] = denoise_dct(params, s, state)
%     params.mask   DCT filter mask
%     state         DSS algorithm state (not used)
%     s             Source signal, matrix of row vector signals
%     s_new         Filtered signal

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'DCT filter';
    params.description = '';
    params.param = {'mask'};
    params.param_value ={[]};
    params.param_type = {'vector'};
    params.param_desc = {'DCT filter coefficients'};
    params.approach = {'pca', 'defl', 'symm'};
    params.beta = {'beta_global'};
    return;
end

% DCT filtering
s_new = idct(repmat(params.mask, size(s,1), 1)' .* dct(s'))';
