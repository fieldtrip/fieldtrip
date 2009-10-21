function [r] = issubfield(s, f)

% ISSUBFIELD tests for the presence of a field in a structure just like the standard
% Matlab ISFIELD function, except that you can also specify nested fields
% using a '.' in the fieldname. The nesting can be arbitrary deep.
%
% Use as
%   f = issubfield(s, 'fieldname')
% or as
%   f = issubfield(s, 'fieldname.subfieldname')
%
% This function returns true if the field is present and false if the field
% is not present.
%
% See also ISFIELD, GETSUBFIELD, SETSUBFIELD

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: issubfield.m,v $
% Revision 1.2  2009/07/30 20:11:44  ingnie
% made output boolian
%
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.1  2005/02/08 09:30:18  roboos
% new implementations to make it easier to work with nested structures
%

try
  getsubfield(s, f);    % if this works, then the subfield must be present  
  r = true;
catch
  r = false;                % apparently the subfield is not present
end