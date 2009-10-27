function [newval, change] = smartinput(question, oldval);

% SMARTINPUT helper function for smart interactive input from the command line
%
% Use as
%   [newval, change] = smartinput(question, oldval)
%
% See also INPUT, PAUSE

% Copyright (C) 2006, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if ischar(oldval)
  newval = input(question, 's');
else
  newval = input(question);
end
if isempty(newval)
  newval = oldval;
  change = 0;
elseif isempty(oldval) && ~isempty(newval)
  change = 1;
elseif ischar(oldval) && strcmp(oldval, newval)
  change = 0;
elseif ~ischar(oldval) && all(size(oldval)==size(newval)) && all(oldval==newval)
  change = 0;
else
  change = 1;
end

