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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

try
  getsubfield(s, f);    % if this works, then the subfield must be present  
  r = true;
catch
  r = false;                % apparently the subfield is not present
end
