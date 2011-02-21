function val = ft_getopt(optarg, key, default)

% FT_GETOPT gets the value of a specified option from a configuration structure
% or from a cell-array with key-value pairs
%
% Use as
%   val = ft_getopt(s, key, default)
% where s is a structure or a cell array. It will return the value of the option, 
% or an empty array if the option was not present. 
%
% See also FT_SETOPT, FT_CHECKOPT

% Copyright (C) 2011, Robert Oostenveld
%
% $Id$

if nargin<3
  default = [];
end

if isa(optarg, 'struct') || isa(optarg, 'config')
  % convert it to a key-value cell-array
  fn = fieldnames(optarg);
  fv = cell(size(fn));
  for i=1:numel(fv)
    fv{i} = optarg.(fn{i});
  end
  optarg = cat(1, fn(:)', fv(:)');
  optarg = optarg(:)';
end

if mod(length(optarg),2)
  error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
end

% the 1st, 3rd, etc. contain the keys, the 2nd, 4th, etc. contain the values
keys = optarg(1:2:end);
vals = optarg(2:2:end);

% the following may be faster than cellfun(@ischar, keys)
valid = false(size(keys));
for i=1:numel(keys)
  valid(i) = ischar(keys{i});
end

if ~all(valid)
  error('optional input arguments should come in key-value pairs, the optional input argument %d is invalid (should be a string)', i);
end

hit = find(strcmpi(key, keys));
if isempty(hit)
  % the requested key was not found
  val = default;
elseif length(hit)==1  
  % the requested key was found
  val = vals{hit};
else
  error('multiple input arguments with the same name');
end
