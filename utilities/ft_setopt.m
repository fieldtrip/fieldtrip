function opt = ft_setopt(opt, key, val)

% FT_SETOPT assigns a value to an configuration structure or to a cell-array
% with key-value pairs. It will overwrite the option if already present, or
% append the option if not present.
%
% Use as
%   s = ft_setopt(s, key, val)
% where s is a structure or a cell array.
%
% See also FT_GETOPT, FT_CHECKOPT

% Copyright (C) 2011-2012, Robert Oostenveld
%
% $Id$

if isa(opt, 'struct') || isa(opt, 'config')
  
  % just replace or add the option
  opt.(key) = val;
  
elseif isa(opt, 'cell')
  % determine whether the key already exists
  fn = opt(1:2:end);
  sel = find(strcmp(key, fn));
  if isempty(sel)
    % append it
    opt{end+1} = key;
    opt{end+1} = val;
  elseif length(sel)==1
    % replace the current value
    keyindex = 2*sel-1;
    valindex = keyindex+1;
    opt{valindex} = val;
  elseif length(sel)>1
    % first remove all occurences
    keyindex = 2*sel-1;
    valindex = keyindex+1;
    opt([keyindex valindex]) = [];
    % then append it
    opt{end+1} = key;
    opt{end+1} = val;
  end
end % isstruct or iscell
