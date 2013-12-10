function dss_message(state, level, message)
% Logging based on verbosity level
%   dss_message(state, level, message)
%     state    DSS state with verbose-parameter
%     level    Used logging level
%     message  Logged message

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if ~isfield(state,'verbose'), return, end
if state.verbose>=level
  fprintf(message);
end
