function M = setdiag(M, v)
% SETDIAG Set the diagonal of a matrix to a specified scalar/vector.
% M = set_diag(M, v)


n = length(M);
if length(v)==1
  v = repmat(v, 1, n);
end
J = 1:n+1:n^2;
M(J) = v;

%M = triu(M,1) + tril(M,-1) + diag(v);

