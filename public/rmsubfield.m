function [s] = rmsubfield(s, f, v);

% RMSUBFIELD removes the contents of the specified field from a structure
% just like the standard Matlab RMFIELD function, except that you can also
% specify nested fields using a '.' in the fieldname. The nesting can be
% arbitrary deep.
%
% Use as
%   s = rmsubfield(s, 'fieldname')
% or as
%   s = rmsubfield(s, 'fieldname.subfieldname')
%
% See also SETFIELD, GETSUBFIELD, ISSUBFIELD

% Copyright (C) 2006, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if ~isstr(f)
  error('incorrect input argument for fieldname');
end

% remove the nested subfield using recursion
[t, f] = strtok(f, '.');
if any(f=='.')
  u = rmsubfield(getfield(s, t), f);
  s = setfield(s, t, u);
else
  s = rmfield(s, t);
end
