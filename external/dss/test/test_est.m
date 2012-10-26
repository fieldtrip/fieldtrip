function test_est
%

% Copyright (C) 2004, 2005 DSS MATLAB package team
% (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of
% Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: test_est.m,v 1.4 2005/05/06 12:43:04 kosti Exp $

T = 2048;
t = 1:T;

S = (1+sin(t*0.01).^4).*randn(1,T);

disp('Gaussian noise:');
mask = estimate_mask(randn(1,T));

sum(mask)

%S = S - mean(S)*ones(1,T);

%input('ENTER','s');

disp('Test signal:');
mask = estimate_mask(S);

sum(mask)
