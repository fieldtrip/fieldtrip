function opt = ft_checkopt(opt, key, allowedtype, allowedval)

% FT_CHECKOPT does a validity test on the types and values of a configuration
% structure or cell-array with key-value pairs.
%
% Use as
%   opt = ft_checkopt(opt, key, allowedtype, allowedval)
%
% For allowedtype you can specify a string or a cell-array with multiple
% strings. All the default MATLAB types can be specified, for example
%   'double'
%   'logical'
%   'char'
%   'single'
%   'float'
%   'int16'
%   'cell'
%   'struct'
%   'function_handle'
% Furthermore, the following custom types can be specified
%   'doublescalar'
%   'doublevector'
%   'doublematrix'
%   'charcell'
%
% For allowedval you can specify a single value or a cell-array
% with multiple values.
%
% This function will give an error or it returns the input configuration
% structure or cell-array without modifications. Any match on allowedtype and
% any match on allowedval is sufficient to let this function pass.
%
% See also FT_GETOPT, FT_SETOPT

% Copyright (C) 2011, Robert Oostenveld
%
% $Id$

if nargin<3
  allowedtype = {};
end

if ~iscell(allowedtype)
  allowedtype = {allowedtype};
end

if nargin<4
  allowedval = {};
end

if ~iscell(allowedval)
  allowedval = {allowedval};
end

% get the value that belongs to this key
val = ft_getopt(opt, key);

% check that the type of the option is allowed
ok = isempty(allowedtype);
for i=1:length(allowedtype)
  switch allowedtype{i}
    case 'empty'
      ok = isempty(val);
    case 'doublescalar'
      ok = isa(val, 'double') && all(size(val)==1);
    case 'doublevector'
      ok = isa(val, 'double') && sum(size(val)>1)==1;
    case 'doublematrix'
      ok = isa(val, 'double') && sum(size(val)>1)>1;
    case 'charcell'
      ok = isa(val, 'cell') && all(cellfun(@ischar, val(:)));
    otherwise
      ok = isa(val, allowedtype{i});
  end
  if ok
    % no reason to do additional checks
    break
  end
end % for allowedtype

if ~ok
  if length(allowedtype)==1
    error('the type of the option "%s" is invalid, it should be "%s" instead of "%s"', key, allowedtype{1}, class(val));
  else
    error('the type of the option "%s" is invalid, it should be any of %s instead of "%s"', key, printcell(allowedtype), class(val));
  end
end

% check that the type of the option is allowed
ok = isempty(allowedval);
for i=1:length(allowedval)
  ok = isequal(val, allowedval{i});
  if ok
    % no reason to do additional checks
    break
  end
end % for allowedtype

if ~ok
  error('the value of the option "%s" is invalid', key);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = printcell(c)
if ~isempty(c)
  s = sprintf('%s, ', c{:});
  s = sprintf('{%s}', s(1:end-2));
else
  s = '{}';
end
