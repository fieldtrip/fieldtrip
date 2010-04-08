function y = any(x)

% ANY    True if any element of a vector is a nonzero number or is
%    logical 1 (TRUE).  ANY ignores entries that are NaN (Not a Number).
%
%    For vectors, ANY(V) returns logical 1 (TRUE) if any of the
%    elements of the vector is a nonzero number or is logical 1 (TRUE).
%    Otherwise it returns logical 0 (FALSE).  For matrices, ANY(X)
%    operates on the columns of X, returning a row vector of logical 1's
%    and 0's.

if nargin>1
  error('this implementation is only supported with one input argument');
end

siz = size(x);
if numel(siz)>2
  error('this implementation is only supported with vector or matrix input');
end

if siz(1)==1 || siz(2)==1
  y = false;
  for i=1:prod(siz)
    if x(i)
      y = true;
      break
    end
  end

else
  y = false(1,siz(2));
  for j=1:siz(2)
    for i=1:siz(1)
      if x(i,j)
        y(1,j) = true;
        break
      end
    end
  end

end



