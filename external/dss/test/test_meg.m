function test_meg

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: test_meg.m,v 1.4 2005/12/07 11:25:17 jaakkos Exp $

% disabled
return;

load '../../meg.mat'

X = xmeg;

T=size(X, 2);
t=1:T;

params.denf.h = @denoise_energy;
params.denf.params.filter_dct = exp(-0.5*(t-1).^2 / (T/64).^2 );

params.reportf = {@test_report_defl @report_w};

params.sdim = 10;

[state, W] = denss(X, params);

S = W * X;

test_report_result(S, state);
