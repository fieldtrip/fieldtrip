function y = diff(x)

% DIFF Difference and approximate derivative.
%    DIFF(X), for a vector X, is [X(2)-X(1)  X(3)-X(2) ... X(n)-X(n-1)].
%    DIFF(X), for a matrix X, is the matrix of row differences,
%       [X(2:n,:) - X(1:n-1,:)].
%    DIFF(X), for an N-D array X, is the difference along the first
%       non-singleton dimension of X.

if nargin>1
  error('this implementation is only supported with one input argument');
end

siz = size(x);
if numel(siz)>2
  error('this implementation is only supported with vector or matrix input');
end

if siz(1)==1
  % derivative along the second dimension
  y = x(:,1:end-1);
  for i=1:(siz(2)-1)
    y(:,i) = x(:,end) - y(:,i);
  end

elseif siz(2)==1
  % derivative along the first dimension
  y = x(1:end-1,:);
  for i=1:(siz(1)-1)
    y(i,:) = x(end,:) - y(i,:);
  end

else
  % derivative along the first dimension
  y = x(1:end-1,:);
  for i=1:(siz(1)-1)
    y(i,:) = x(end,:) - y(i,:);
  end
end

