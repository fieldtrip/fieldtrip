function [report_data] = report_convergence(report_data, state)
% Report convergence of the algorithm

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if ~isfield(report_data, 'report_interval')
  report_data.report_interval = 5;
  dss_message(state, 2, sprintf('Setting default report interval to %d (report_data.report_interval)\n', report_data.report_interval));
end

% deflation
if state.iteration==1
  % -- first iteration
  if state.algorithm=='defl'
    % deflation
    report_data.change(state.component,1)=0;
    report_data.deltaw_old=state.w;
  else
    % symmetric
    report_data.change = zeros(size(state.W, 1),1);
    report_data.dW_old = zeros(size(state.W));
    report_data.W_old = state.W;
  end
else
  % -- iterations 2->
  if (mod(state.iteration, report_data.report_interval)==0)
    if state.algorithm=='defl'
      % deflation
      change = acos(state.w' * state.w_old/norm(state.w)/norm(state.w_old)) / pi * 180;
    else 
      % symmetric
      change = abs(angle(state.W, state.W_old)) / pi * 180;
    end
    message(state,1,sprintf('Change (angles): %d\n', change));
  end
end

% -----------------------------------
function rad = angle(A, B)

sum_cross = sum(A .* B, 2);
sum_A = sum(A .* A, 2);
sum_B = sum(B .* B, 2);
rad = acos(sum_cross ./ sum_A.^(-1/2) ./ sum_B.^(-1/2));
