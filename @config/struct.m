function y = struct(x, varargin)

% STRUCT Convert a config object into a structure object.

if nargin==1
  % convert the config object into a regular Matlab structure
  for i=1:numel(x)
    y(i) = struct(x(i).value);
  end
  y = reshape(y, size(x));
  % recurse into the structure to convert sub-configs into sub-structures
  key = fieldnames(y);
  for i=1:length(key)
    val = y.(key{i});
    if isa(val, 'config')
      y = setfield(y, key{i}, struct(val));
    end
  end
else
  % mimic the behaviour of the builtin Matlab struct function
  if mod(nargin,2)
    error('Incorrect number of input arguments (should be key-value pairs)')
  end
  varargin = {x varargin{:}};
  key = varargin(1:2:end);
  val = varargin(2:2:end);

  y = struct();
  for i=1:length(key)
    y = setfield(y, key{i}, val{i});;
  end
end
