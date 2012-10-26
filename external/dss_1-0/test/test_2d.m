function test_2d

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: test_2d.m,v 1.5 2005/12/07 11:25:17 jaakkos Exp $

X = create_mixed_source(4, 256*256);

X = reshape(X, 4, 256, []);

state = denss(X);

% convert result components back to 2 dimensions
S = reshape(state.S, [size(state.S,1), state.input_dims]);

s = size(S);
if (length(s)~=3) | any(s~=[4 256 256])
    error('Invalid result component dimension');
end