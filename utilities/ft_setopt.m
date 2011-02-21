function opt = ft_setopt(opt, key, val)

% FT_SETOPT assigns a value to an configuration structure or to a
% cell-array with key-value pairs. It will overwrite the option if already
% present, or append the option if not present.
%
% Use as
%   opt = ft_setopt(s, key, val)
% where s is a structure or a cell array.
%
% See also T_GETOPT, FT_CHECKOPT

% Copyright (C) 2011, Robert Oostenveld
%
% $Id$

if isa(optarg, 'struct') || isa(optarg, 'config')
  % just replace or add the option
  opt.(key) = val;
elseif iscell
  fn = opt(1:2:end);
  sel = find(strcmp(key, fn));
  if isempty(sel)
    % append it
    opt{end+1} = key;
    opt{end+1} = val;
  else
    % replace the current value
    opt{sel+1} = val;
  end
end
