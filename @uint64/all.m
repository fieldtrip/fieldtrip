function y = any(x)

% ALL    True if all elements of a vector are nonzero.
%    For vectors, ALL(V) returns logical 1 (TRUE) if none of the elements 
%    of the vector are zero.  Otherwise it returns logical 0 (FALSE).  For 
%    matrices, ALL(X) operates on the columns of X, returning a row vector
%    of logical 1's and 0's. For N-D arrays, ALL(X) operates on the first
%    non-singleton dimension.

if nargin>1
  error('this implementation is only supported with one input argument');
end

siz = size(x);
if numel(siz)>2
  error('this implementation is only supported with vector or matrix input');
end

if siz(1)==1 || siz(2)==1
  y = true;
  for i=1:prod(siz)
    if ~x(i)
      y = false;
      break
    end
  end

else
  y = true(1,siz(2));
  for j=1:siz(2)
    for i=1:siz(1)
      if ~x(i,j)
        y(1,j) = false;
        break
      end
    end
  end

end



