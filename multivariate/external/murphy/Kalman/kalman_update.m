function [xnew, Vnew, loglik, VVnew] = kalman_update(A, C, Q, R, y, x, V, varargin)
% KALMAN_UPDATE Do a one step update of the Kalman filter
% [xnew, Vnew, loglik] = kalman_update(A, C, Q, R, y, x, V, ...)
%
% INPUTS:
% A - the system matrix
% C - the observation matrix 
% Q - the system covariance 
% R - the observation covariance
% y(:)   - the observation at time t
% x(:) - E[X | y(:, 1:t-1)] prior mean
% V(:,:) - Cov[X | y(:, 1:t-1)] prior covariance
%
% OPTIONAL INPUTS (string/value pairs [default in brackets])
% 'initial' - 1 means x and V are taken as initial conditions (so A and Q are ignored) [0]
% 'u'     - u(:) the control signal at time t [ [] ]
% 'B'     - the input regression matrix
%
% OUTPUTS (where X is the hidden state being estimated)
%  xnew(:) =   E[ X | y(:, 1:t) ] 
%  Vnew(:,:) = Var[ X(t) | y(:, 1:t) ]
%  VVnew(:,:) = Cov[ X(t), X(t-1) | y(:, 1:t) ]
%  loglik = log P(y(:,t) | y(:,1:t-1)) log-likelihood of innovatio

% set default params
u = [];
B = [];
initial = 0;

args = varargin;
for i=1:2:length(args)
  switch args{i}
   case 'u', u = args{i+1};
   case 'B', B = args{i+1};
   case 'initial', initial = args{i+1};
   otherwise, error(['unrecognized argument ' args{i}])
  end
end

%  xpred(:) = E[X_t+1 | y(:, 1:t)]
%  Vpred(:,:) = Cov[X_t+1 | y(:, 1:t)]

if initial
  if isempty(u)
    xpred = x;
  else
    xpred = x + B*u;
  end
  Vpred = V;
else
  if isempty(u)
    xpred = A*x;
  else
    xpred = A*x + B*u;
  end
  Vpred = A*V*A' + Q;
end

e = y - C*xpred; % error (innovation)
n = length(e);
ss = length(A);
S = C*Vpred*C' + R;
Sinv = inv(S);
ss = length(V);
loglik = gaussian_prob(e, zeros(1,length(e)), S, 1);
K = Vpred*C'*Sinv; % Kalman gain matrix
% If there is no observation vector, set K = zeros(ss).
xnew = xpred + K*e;
Vnew = (eye(ss) - K*C)*Vpred;
VVnew = (eye(ss) - K*C)*A*V;

