function [coef, C] = learn_AR(data, k)
% Find the ML parameters of a vector autoregressive process of order k.
% [coef, C] = learn_AR(k, data)
% data{l}(:,t) = the observations at time t in sequence l

warning('learn_AR seems to be broken');

nex = length(data);
obs = cell(1, nex);
for l=1:nex
  obs{l} = convert_to_lagged_form(data{l}, k);
end

% The initial parameter values don't matter, since this is a perfectly observable problem.
% However, the size of F must be set correctly.
y = data{1};
[s T] = size(y);
coef = rand(s,s,k);
C = rand_psd(s);
[F,H,Q,R,initx,initV] = AR_to_SS(coef, C, y);

max_iter = 1;
fully_observed = 1;
diagQ = 0;
diagR = 0;
[F, H, Q, R, initx, initV, loglik] = ...
    learn_kalman(obs, F, H, Q, R, initx, initV, max_iter, diagQ, diagR, fully_observed);

[coef, C] = SS_to_AR(F, Q, k);

