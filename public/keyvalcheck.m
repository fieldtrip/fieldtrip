function keyvalcheck(arglist, varargin)

% KEYVALCHECK is a helper function for parsing optional key-value input pairs.
%
% Use as
%   keyvalcheck(argin, 'required',  {'key1', 'key2', ...})
%   keyvalcheck(argin, 'forbidden', {'key1', 'key2', ...})
%   keyvalcheck(argin, 'optional',  {'key1', 'key2', ...})
%
% See also KEYVAL

% Copyright (C) 2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% this is to speed up subsequent calls with the same input arguments
persistent previous_argin

current_argin = {arglist, varargin{:}};
if ~isempty(previous_argin) && isequal(previous_argin, current_argin)
  % the input is the same to the previous input, and that was OK
  return
end

required  = keyval('required', varargin);
forbidden = keyval('forbidden', varargin);
optional  = keyval('optional', varargin);

keys = arglist(1:2:end);
vals = arglist(2:2:end);

if numel(keys)~=numel(vals)
  error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
end

keys = cellfun(@lower, keys, 'UniformOutput', false);

if ~isempty(required)
  % only check if specified
  required  = cellfun(@lower, required, 'UniformOutput', false);
  set = intersect(keys, required);
  if numel(set)~=numel(required)
    error('the required input argument ''%s'' was not specified', set{:});
  end
end

if ~isempty(forbidden)
  % only check if specified
  forbidden = cellfun(@lower, forbidden, 'UniformOutput', false);
  set = intersect(keys, forbidden);
  if numel(set)~=0
    error('the input argument ''%s'' is forbidden', set{:});
  end
end

if ~isempty(optional)
  % only check if specified
  optional  = cellfun(@lower, optional, 'UniformOutput', false);
  set = setdiff(keys, optional);
  if numel(set)>0
    error('the input argument ''%s'' is forbidden', set{:});
  end
end

% remember the current input arguments, which appear to be OK
previous_argin = current_argin;
