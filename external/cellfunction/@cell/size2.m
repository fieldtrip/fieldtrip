function [varargout] = size2(x, dim, cellflag)

% SIZE2 Size of cells in cell-array.
% 
% Use as:
%   [varargout] = size2(x, dim, cellflag)
%
% Input arguments:
%   x   = cell-array
%   dim = dimension argument, default []
%   cellflag = true (default) or false. If false, size2 behaves as size,
%     otherwise the size of each of the cells is returned.

if nargin<3, cellflag = false; end
if nargin<2, dim      = [];    end

if nargout>1 && ~isempty(dim)
  error('Too many output arguments.');
end

if ischar(cellflag) && strcmp(cellflag, 'cell')
  cellflag = true;
end

if ~cellflag
  % traditional functionality
  numdim = ndims(x);
  
  if isempty(dim)
    siz    = zeros(1,numdim);
    for k = 1:numdim
      switch k
        case 1
          siz(k) = numel(x(:,1,1,1,1,1,1,1,1,1));
        case 2
          siz(k) = numel(x(1,:,1,1,1,1,1,1,1,1));
        case 3
          siz(k) = numel(x(1,1,:,1,1,1,1,1,1,1));
        case 4
          siz(k) = numel(x(1,1,1,:,1,1,1,1,1,1));
        case 5
          siz(k) = numel(x(1,1,1,1,:,1,1,1,1,1));
        case 6
          siz(k) = numel(x(1,1,1,1,1,:,1,1,1,1));
        case 7
          siz(k) = numel(x(1,1,1,1,1,1,:,1,1,1));
        case 8
          siz(k) = numel(x(1,1,1,1,1,1,1,:,1,1));
        case 9
          siz(k) = numel(x(1,1,1,1,1,1,1,1,:,1));
        case 10
          siz(k) = numel(x(1,1,1,1,1,1,1,1,1,:));
        otherwise
      end
    end
  else
    switch dim
      case 1
        siz = numel(x(:,1,1,1,1,1,1,1,1,1));
      case 2
        siz = numel(x(1,:,1,1,1,1,1,1,1,1));
      case 3
        siz = numel(x(1,1,:,1,1,1,1,1,1,1));
      case 4
        siz = numel(x(1,1,1,:,1,1,1,1,1,1));
      case 5
        siz = numel(x(1,1,1,1,:,1,1,1,1,1));
      case 6
        siz = numel(x(1,1,1,1,1,:,1,1,1,1));
      case 7
        siz = numel(x(1,1,1,1,1,1,:,1,1,1));
      case 8
        siz = numel(x(1,1,1,1,1,1,1,:,1,1));
      case 9
        siz = numel(x(1,1,1,1,1,1,1,1,:,1));
      case 10
        siz = numel(x(1,1,1,1,1,1,1,1,1,:));
      otherwise
    end
  end
  
  if nargout<=1
    varargout{1} = siz;
  elseif nargout>=numdim
    for k = 1:numdim,           varargout{k} = siz(k); end
    for k = (numdim+1):nargout, varargout{k} = 1;      end
  elseif nargout<numdim
    for k = 1:(nargout-1),      varargout{k} = siz(k); end
    varargout{nargout} = prod(siz((k+1):end));
  end
  return;
else
  % operate on the individual cells
  siz    = ones(numel(x), 10);
  numdim = zeros(numel(x), 1);
  for k = 1:numel(x)
    numdim(k) = ndims(x{k});
    siz(k,1:numdim(k)) = size(x{k});
  end
 
  if isempty(dim) 
    if nargout<=1
      varargout{1} = siz(:,1:max(numdim));
    elseif nargout>=max(numdim)
      for k = 1:numdim,           varargout{k} = reshape(siz(:,k), size(x)); end
      for k = (numdim+1):nargout, varargout{k} = ones(size(x));              end  
    end
  else
    varargout{1} = reshape(siz(:,dim), size(x));
  end
  return;
end

  


