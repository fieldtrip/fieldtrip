function X = create_mixed_source(mixed_dim, T, seed)

% Create pseudorandom mixed signals

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: create_mixed_source.m,v 1.4 2005/05/24 08:42:40 kosti Exp $

% Number of mixed signals
M = mixed_dim;

% Signal lengh
if nargin<2; T=2048; end;
if nargin<3; seed=0; end;
t=1:T;

randn('state',seed);

S = [sin(t*0.02) ; sin(t*0.07) ; sin(t*0.27) ; randn(1,T) ; sin(t*0.017).^9 ; mod(t, 64)/6];

% Number of original signals
N = size(S,1);

A = randn(M,N);
X = A * S;
