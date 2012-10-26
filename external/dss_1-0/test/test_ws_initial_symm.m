function test_ws_initial_symm
% Test initial guesses for W and S with symmetric DSS

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: test_ws_initial_symm.m,v 1.2 2005/12/07 11:25:17 jaakkos Exp $

X = create_mixed_source(8);

[dummy, Y] = pre_sphere([], X);

T = size(Y,2);
M = size(Y,1);

params.Y = Y;
params.algorithm = 'symm';
params.W = randn(2, M);
params.S = randn(4, T);

state = dss_create_state(params);
%diary('temp_output')
state = denss(state);
%diary off
