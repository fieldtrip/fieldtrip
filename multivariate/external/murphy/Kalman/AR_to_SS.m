function [F,H,Q,R,initx, initV] = AR_to_SS(coef, C, y)
%
% Convert a vector auto-regressive model of order k to state-space form.
% [F,H,Q,R] = AR_to_SS(coef, C, y)
% 
% X(i) = A(1) X(i-1) + ... + A(k) X(i-k+1) + v, where v ~ N(0, C)
% and A(i) = coef(:,:,i) is the weight matrix for i steps ago.
% We initialize the state vector with [y(:,k)' ... y(:,1)']', since
% the state vector stores [X(i) ... X(i-k+1)]' in order.

[s s2 k] = size(coef); % s is the size of the state vector
bs = s * ones(1,k); % size of each block

F = zeros(s*k);
for i=1:k
   F(block(1,bs), block(i,bs)) = coef(:,:,i);
end
for i=1:k-1
  F(block(i+1,bs), block(i,bs)) = eye(s);
end

H = zeros(1*s, k*s);
% we get to see the most recent component of the state vector 
H(block(1,bs), block(1,bs)) = eye(s); 
%for i=1:k
%  H(block(1,bs), block(i,bs)) = eye(s);
%end

Q = zeros(k*s);
Q(block(1,bs), block(1,bs)) = C;

R = zeros(s);

initx = zeros(k*s, 1);
for i=1:k
  initx(block(i,bs)) = y(:, k-i+1); % concatenate the first k observation vectors
end

initV = zeros(k*s); % no uncertainty about the state (since perfectly observable)
