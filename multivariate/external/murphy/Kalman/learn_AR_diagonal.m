function [coef, C] = learn_AR_diagonal(y, k)
% Find the ML parameters for a collection of independent scalar AR processes.

% sep_coef(1,1,t,i) is the coefficient to apply to compopnent i of the state vector t steps ago
% eg. consider two components L and R and let A = coef(:,:,1,:), B = coef(:,:,2,:)
% L3    (AL  0   BL 0)     (L2)   (CL 0 0 0)
% R3 =  (0   AR  0  BR)    (R2)   (0 CR 0 0)
% L2    (1   0   0  0 )    (L1) + (0 0  0 0)
% R2    (0   1   0  0 )    (R1)   (0 0  0 0)

ss = size(y, 1);
sep_coef = zeros(1,1,k,ss);
for i=1:ss
  [sep_coef(:,:,:,i), sep_cov(i)] = learn_AR(k, y(i,:));
end
C = diag(sep_cov);
for t=1:k
  x = sep_coef(1,1,t,:);
  coef(:,:,t) = diag(x(:));
end
