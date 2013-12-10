function [report_data] = report_print(report_data, state)
% Example of iteration reporting function.
% Prints some information about the component.
%

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

data = state;

if data.iteration==1
  % -- first iteration
  dss_message(state,1,sprintf('  Iterating component %d', data.component));
end

% -- show progress status to the user
dss_message(state,1,'.');
if data.last_iteration
  dss_message(state,1,sprintf('\n  Component %d found after %d iterations.\n', ...
	  data.component, data.iteration));

  dss_message(state,1,sprintf('  mean: %f, var: %f, skew: %f, kurt: %f\n', ...
	  mean(data.s), var(data.s,1), skewness(data.s), kurtosis(data.s)-3));
end

