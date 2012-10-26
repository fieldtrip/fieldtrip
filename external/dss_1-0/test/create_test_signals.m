function X = create_test_signals(signal_length)
% Create artificial signals for testing.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: create_test_signals.m,v 1.2 2005/04/20 10:19:25 kosti Exp $

% Number of original signals
N=6;

% Number of mixed signals
M=10;

% Signal lentgh
T=signal_length;
t=1:T;

% Create random signals
S_orig=[sin(t*0.02) ; randn(1,T) ; randlap(1,T) ; sin(t*0.017).^9 ; mod(t, 64)/6 ; sin((t-T/2).^2*0.0005)];

% -- Normalize S
S_orig = S_orig-repmat(mean(S_orig')',1,T);
S_orig = diag(std(S_orig', 1).^(-1))*S_orig;

% -- Create mixed signals
A=randn(M,N);
X=A*S_orig+randn(M,T)*0.5;
