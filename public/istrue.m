function y = istrue(x)

% ISTRUE ensures that a true/false input argument like "yes", "true"
% or "on" is converted into a boolean

% Copyright (C) 2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

