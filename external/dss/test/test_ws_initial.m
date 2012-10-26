function test_ws_initial
% Test initial guesses for W and S
%
% Test should output following lines:
%   Calculate component 1 with predefined w
%   Calculate component 2 with predefined w
%   Calculate component 3 with predefined s
%   Calculate component 4 with predefined s
%   Calculate component 5 with random w
%   Calculate component 6 with random w

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: test_ws_initial.m,v 1.2 2005/12/07 11:25:17 jaakkos Exp $

X = create_mixed_source(8);

[dummy, Y] = pre_sphere([], X);

T = size(Y,2);
M = size(Y,1);

params.Y = Y;
params.algorithm = 'defl';
params.W = ones(2, M);
params.S = ones(4, T);

state = dss_create_state(params);
%diary('temp_output')
state = denss(state);
%diary off
