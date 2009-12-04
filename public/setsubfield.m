function [s] = setsubfield(s, f, v);

% SETSUBFIELD sets the contents of the specified field to a specified value
% just like the standard Matlab SETFIELD function, except that you can also
% specify nested fields using a '.' in the fieldname. The nesting can be
% arbitrary deep.
%
% Use as
%   s = setsubfield(s, 'fieldname', value)
% or as
%   s = setsubfield(s, 'fieldname.subfieldname', value)
%
% See also SETFIELD, GETSUBFIELD, ISSUBFIELD

% Copyright (C) 2005, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if ~isstr(f)
  error('incorrect input argument for fieldname');
end

t = {};
while (1)
  [t{end+1}, f] = strtok(f, '.');
  if isempty(f)
    break
  end
end
s = setfield(s, t{:}, v);
