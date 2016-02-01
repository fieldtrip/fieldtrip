% Prepares DSS calculation environment. Currently only compiles the
% testkeypress MEX function.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

disp('Initializing DSS environment');
disp('Compiling testkeypress executable');

mex src/testkeypress.c;

disp('Environment initialized successfully');
