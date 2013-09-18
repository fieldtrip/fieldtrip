function [errorcode, varargout] = common_size (varargin)

if (nargin < 2)
  error ('common_size: only makes sense if nargin >= 2');
end

% Find scalar args.
nscal = cellfun (@numel, varargin) ~= 1;

i = find (nscal, 1);

if (isempty (i))
  errorcode = 0;
  varargout = varargin;
else
  match = cellfun (@size_equal, varargin, repmat(varargin(i),size(varargin)));
  if (any (nscal & ~match))
    errorcode = 1;
    varargout = varargin;
  else
    errorcode = 0;
    if (nargout > 1)
      scal = ~nscal;
      varargout = varargin;
      dims = size (varargin{i});
      for s = find(scal)
        varargout{s} = repmat(varargin{s}, dims);
      end
    end
  end
end
end