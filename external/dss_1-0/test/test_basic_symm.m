function test_basic_symm

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: test_basic_symm.m,v 1.3 2005/12/07 11:25:17 jaakkos Exp $

X = create_mixed_source(4);

params.algorithm = 'symm';

state = denss(X, params);

