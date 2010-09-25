function [coef, C] = SS_to_AR(F, Q, k, diagonal)
%
% Extract the parameters of a vector autoregresssive process of order k from the state-space form.
% [coef, C] = SS_to_AR(F, Q, k, diagonal)

if nargin<4, diagonal = 0; end

s = length(Q) / k;
bs = s*ones(1,k);
coef = zeros(s,s,k);
for i=1:k
  if diagonal
    coef(:,:,i) = diag(diag(F(block(1,bs), block(i,bs))));
  else
    coef(:,:,i) = F(block(1,bs), block(i,bs));
  end
end
C = Q(block(1,bs), block(1,bs));
if diagonal
  C = diag(diag(C));
end
%C = sqrt(Q(block(1,bs), block(1,bs))); % since cov(1,1) of full vector = C C'
