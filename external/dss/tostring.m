function [s] = tostring(values)
% convert parameter value to string

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if strcmp(class(values), 'char')
  s = values;
  return;
end
if size(values,1)>1 & size(values,2)>1
    s = 'matrix';
    return;
end
s = '';
for v=1:size(values,2)
  value = values(v);
  if strcmp(class(value),'cell');
      value = value{1};
  end
  if v>1; s = [s ', ']; end
  switch class(value)
   case 'struct'
    if isfield(value, 'h')
      s = [s '@' func2str(value.h)];
    else
      s = [s 'structure'];
    end
   case 'logical'
    if value s = [s 'true'];
    else s = [s 'false'];
    end;
   case 'function_handle'
    s = [s '@' func2str(value)];
   case 'double'
    s = [s num2str(value)];
   otherwise
    s = [s value];
  end
end
  