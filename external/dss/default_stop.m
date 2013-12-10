function [stop, params] = default_stop(params, state)
% Default stopping criteria function. 
%   [stop, params] = default_stop(params, state)
%     params  Function specific modifiable parameters
%     params.maxiters  Maximum number of allowed iterations
%     params.epsilon   Treshold angle, iteration is stopped when
%                      angle beween consecutive iteration
%                      projection vectors goes below treshold.
%     state   DSS algorithm state
%     stop    Boolean, true when stopping criteria has been met

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if ~isfield(params, 'maxiters') & ...
      (~isfield(params, 'epsilon') | ~isfield(state, 'w'))
  error('Either ''maxiters'' or ''epsilon'' must be defined as stopping criteria');
end

stop = 0;

if isfield(params, 'maxiters')
  if state.iteration>=params.maxiters
    stop = 1;
    return;
  end
end

if isfield(params, 'epsilon')
  if isfield(state, 'w')
    change = angle(state.w_old, state.w) / pi * 180;
  else
    change = angle(state.W_old, state.W) / pi * 180;
  end
  stop = all(change < params.epsilon);
end

% --------
function rad = angle(A, B)

sum_cross = sum(A .* B);
sum_A = sum(A .* A);
sum_B = sum(B .* B);
rad = acos(sum_cross .* sum_A.^(-1/2) .* sum_B.^(-1/2));
