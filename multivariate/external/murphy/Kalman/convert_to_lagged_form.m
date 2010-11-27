function yy = convert_to_lagged_form(y, k)
% Create an observation vector yy(:,t) containing the last k values of y, newest first
% e.g., k=2, y = (a1 a2 a3)     yy  = a2 a3
%                (b1 b2 b3)           b2 b2
%                                     a1 a2
%                                     b1 b2

[s T] = size(y);
bs = s*ones(1,k);
yy = zeros(k*s, T-k+1);
for i=1:k
  yy(block(i,bs), :) = y(:, k-i+1:end-i+1);
end

