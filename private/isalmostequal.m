function [ok, message] = isalmostequal(a, b, varargin)

% ISALMOSTEQUAL compares two input variables and returns true/false
% and a message containing the details on the observed difference.
%
% Use as
%   [ok, message] = isalmostequal(a, b)
%   [ok, message] = isalmostequal(a, b, ...)
%
% This works for all possible input variables a and b, like
% numerical arrays, string arrays, cell arrays, structures
% and nested data types.
%
% Optional input arguments come in key-value pairs, supported are
%   'depth'      number, for nested structures
%   'abstol'     number, absolute tolerance for numerical comparison
%   'reltol'     number, relative tolerance for numerical comparison
%   'diffabs'    boolean, check difference between absolute values for numericals (useful for e.g. mixing matrices which have arbitrary signs)
%
% See also ISEQUAL, ISEQUALNAN

% Copyright (C) 2004-2012, Robert Oostenveld & Markus Siegel
%
% $Id$

if nargin==3
  % for backward compatibility
  depth = varargin{1};
else
  depth = ft_getopt(varargin, 'depth');
  if isempty(depth)
    % set the default
    depth = inf;
  end
end

message = {};
location = '';

[message] = do_work(a, b, depth, location, message, varargin{:});
message = message(:);
ok = isempty(message);

if ~nargout
  disp(message);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [message] = do_work(a, b, depth, location, message, varargin)

knowntypes = {
  'double'          % Double precision floating point numeric array
  'logical'         % Logical array
  'char'            % Character array
  'cell'            % Cell array
  'struct'          % Structure array
  'numeric'         % Integer or floating-point array
  'single'          % Single precision floating-point numeric array
  'int8'            % 8-bit signed integer array
  'uint8'           % 8-bit unsigned integer array
  'int16'           % 16-bit signed integer array
  'uint16'          % 16-bit unsigned integer array
  'int32'           % 32-bit signed integer array
  'uint32'          % 32-bit unsigned integer array
  };

for type=knowntypes(:)'
  if isa(a, type{:}) && ~isa(b, type{:})
    message{end+1} = sprintf('different data type in %s', location);
    return
  end
end

if isempty(location)
  location = 'array';
end

if isa(a, 'numeric') || isa(a, 'char') || isa(a, 'logical')
  % perform numerical comparison
  if length(size(a))~=length(size(b))
    message{end+1} = sprintf('different number of dimensions in %s', location);
    return;
  end
  if any(size(a)~=size(b))
    message{end+1} = sprintf('different size in %s', location);
    return;
  end
  if ~all(isnan(a(:)) == isnan(b(:)))
    message{end+1} = sprintf('different occurence of NaNs in %s', location);
    return;
  end
  % replace the NaNs, since we cannot compare them numerically
  a = a(~isnan(a(:)));
  b = b(~isnan(b(:)));
  % continue with numerical comparison
  if ischar(a) && any(a~=b)
    message{end+1} = sprintf('different string in %s: %s ~= %s', location, a, b);
  else
    % use the desired tolerance
    reltol     = ft_getopt(varargin, 'reltol');       % any value, relative to the mean
    abstol     = ft_getopt(varargin, 'abstol');       % any value
    relnormtol = ft_getopt(varargin, 'relnormtol');   % the matrix norm, relative to the mean norm
    absnormtol = ft_getopt(varargin, 'absnormtol');   % the matrix norm
    diffabs    = ft_getopt(varargin, 'diffabs');
    
    if ~isempty(diffabs) && diffabs
      a = abs(a);
      b = abs(b);
    end
    
    if ~isempty(abstol) && any(abs(a-b)>abstol)
      message{end+1} = sprintf('different values in %s', location);
    elseif ~isempty(reltol) && any((abs(a-b)./(0.5*abs(a+b)))>reltol)
      message{end+1} = sprintf('different values in %s', location);
    elseif isempty(abstol) && isempty(reltol) && any(a~=b)
      message{end+1} = sprintf('different values in %s', location);
    elseif ~isempty(relnormtol) && (norm(a-b)/(0.5*(norm(a)+norm(b)))>relnormtol)
      message{end+1} = sprintf('different values in %s', location);
    elseif ~isempty(absnormtol) && norm(a-b)>absnormtol
      message{end+1} = sprintf('different values in %s', location);
    end
  end
  
elseif isa(a, 'struct') && all(size(a)==1)
  % perform recursive comparison of all fields of the structure
  fna = fieldnames(a);
  fnb = fieldnames(b);
  if ~all(ismember(fna, fnb))
    tmp = fna(~ismember(fna, fnb));
    for i=1:length(tmp)
      message{end+1} = sprintf('field missing in the 2nd argument in %s: {%s}', location, tmp{i});
    end
  end
  if ~all(ismember(fnb, fna))
    tmp = fnb(~ismember(fnb, fna));
    for i=1:length(tmp)
      message{end+1} = sprintf('field missing in the 1st argument in %s: {%s}', location, tmp{i});
    end
  end
  fna = intersect(fna, fnb);
  if depth>0
    % warning, this is a recursive call to transverse nested structures
    for i=1:length(fna)
      fn = fna{i};
      ra = getfield(a, fn);
      rb = getfield(b, fn);
      [message] = do_work(ra, rb, depth-1, [location '.' fn], message, varargin{:});
    end
  end
  
elseif isa(a, 'struct') && ~all(size(a)==1)
  % perform recursive comparison of all array elements
  if any(size(a)~=size(b))
    message{end+1} = sprintf('different size of struct-array in %s', location);
    return;
  end
  siz = size(a);
  dim = ndims(a);
  a = a(:);
  b = b(:);
  for i=1:length(a)
    ra = a(i);
    rb = b(i);
    tmp = sprintf('%s(%s)', location, my_ind2sub(siz, i));
    [message] = do_work(ra, rb, depth-1, tmp, message, varargin{:});
  end
  
elseif isa(a, 'cell')
  % perform recursive comparison of all array elements
  if any(size(a)~=size(b))
    message{end+1} = sprintf('different size of cell-array in %s', location);
    return;
  end
  siz = size(a);
  dim = ndims(a);
  a = a(:);
  b = b(:);
  for i=1:length(a)
    ra = a{i};
    rb = b{i};
    tmp = sprintf('%s{%s}', location, my_ind2sub(siz, i));
    [message] = do_work(ra, rb, depth-1, tmp, message, varargin{:});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return a string with the formatted subscript
function [str] = my_ind2sub(siz,ndx)
n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1,
  tmp(i) = floor(ndx/k(i))+1;
  ndx = rem(ndx,k(i));
end
str = '';
for i=1:n
  str = [str ',' num2str(tmp(i))];
end
str = str(2:end);
