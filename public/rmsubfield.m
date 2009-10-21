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
% $Log: rmsubfield.m,v $
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.1  2006/11/27 15:38:45  roboos
% new implementation
%

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
