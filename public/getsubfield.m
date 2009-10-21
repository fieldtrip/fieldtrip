function [s] = getsubfield(s, f);

% GETSUBFIELD returns a field from a structure just like the standard
% Matlab GETFIELD function, except that you can also specify nested fields
% using a '.' in the fieldname. The nesting can be arbitrary deep.
%
% Use as
%   f = getsubfield(s, 'fieldname')
% or as
%   f = getsubfield(s, 'fieldname.subfieldname')
%
% See also GETFIELD, ISSUBFIELD, SETSUBFIELD

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: getsubfield.m,v $
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.1  2005/02/08 09:30:18  roboos
% new implementations to make it easier to work with nested structures
%

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
s = getfield(s, t{:});
