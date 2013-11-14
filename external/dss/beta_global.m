function [params, beta] = beta_global(params, state, denf)
% Beta function for global spectral shift based on gaussian noise.
%   [s_new, params] = beta_global(params, state, s)
%     params  Function specific modifiable parameters
%     state   DSS algorithm state
%     denf    Denoising function
%     beta    Global value for denoising beta

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$
                                                        
if nargin<2
    params.name = 'Global spectral shift for all denoising functions';
    params.adaptive = 0;
    params.approach = {'defl','symm'};
    return;
end

signal_length = size(state.Y, 2);
v = randn(1,signal_length);
beta = -feval(denf.h, denf.params, state, v) * v' /signal_length;
