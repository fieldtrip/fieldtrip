function [params, s_new] = denoise_mask(params, s, state)
%  Mask DSS denoising function
%   [params, s_new] = denoise_mask(params, s, state)
%     params  Function specific modifiable parameters
%     params.mask ...
%     state   DSS algorithm state
%     s       Source signal estimate, matrix of row vector signals
%     s_new   Denoised signal estimate
%
%   Multiplies the signal with simple mask. Can be used for ON/OFF
%   denoising with 0 and 1 values.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'Mask denoising';
    params.description = '';
    params.param = {'mask'};
    params.param_value ={[]};
    params.param_type = {'vector'};
    params.param_desc = {'Denoising mask'};
    params.approach = {'pca','defl','symm'};
    params.alpha = {};
    params.beta = {'beta_global'};
    return;
end

% masking in time
if strcmp(state.algorithm,'defl') 
  s_new = s .* params.mask;
else
  [sdim,tdim] = size(s);
  s_new = s .* repmat(params.mask,sdim,1);
end
