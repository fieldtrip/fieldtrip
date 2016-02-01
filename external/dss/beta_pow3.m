function [params, beta] = beta_pow3(params, state, denf)
% Local spectral shift function for denoise_tanh denoising function.
%   [s_new, params] = beta_pow3(params, state, s)
%     params  Function specific modifiable parameters
%     state   DSS algorithm state
%     denf    Tanh denoising function
%     beta    Denoising beta value

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<2
    params.name = 'Local spectral shift for pow3';
    params.adaptive=1;
    params.approach = {'defl','symm'};
    return;
end

beta = -3;
