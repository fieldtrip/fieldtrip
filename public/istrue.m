function y = istrue(x)

% ISTRUE ensures that a true/false input argument like "yes", "true"
% or "on" is converted into a boolean

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: istrue.m,v $
% Revision 1.4  2009/06/05 14:26:23  crimic
% added 'none' option
%
% Revision 1.3  2009/05/14 11:57:37  crimic
% introduced check for 'y' and 'n'
%
% Revision 1.2  2009/04/14 18:30:35  roboos
% small fix
%
% Revision 1.1  2009/04/14 18:28:45  roboos
% extended and moved from plotting to public
%

true_list  = {'yes' 'true' 'on' 'y' };
false_list = {'no' 'false' 'off' 'n' 'none'};

if ischar(x)
  % convert string to boolean value
  if any(strcmpi(x, true_list))
    y = true;
  elseif any(strcmpi(x, false_list))
    y = false;
  else
    error('cannot determine whether "%s" should be interpreted as true or false', x);
  end
else
  % convert numerical value to boolean
  y = logical(x);
end

