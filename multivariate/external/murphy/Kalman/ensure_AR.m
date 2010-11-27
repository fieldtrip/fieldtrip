function [A, C, Q, R, initx, initV] = ensure_AR(A, C, Q, R, initx, initV, k, obs, diagonal)
%
% Ensure that the system matrices have the right form for an autoregressive process.

ss = length(A);
if nargin<8, obs=ones(ss, 1); end
if nargin<9, diagonal=0; end

[coef, C] = SS_to_AR(A, Q, k, diagonal);
[A, C, Q, R, initx, initV] = AR_to_SS(coef, C, obs);
