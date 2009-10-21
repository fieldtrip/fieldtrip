function [newval, change] = smartinput(question, oldval);

% SMARTINPUT helper function for smart interactive input from the command line
%
% Use as
%   [newval, change] = smartinput(question, oldval)
%
% See also INPUT, PAUSE

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: smartinput.m,v $
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.4  2006/06/01 12:52:21  roboos
% allow for vector and matrix input
%
% Revision 1.3  2006/05/02 19:13:14  roboos
% allow empty oldvalue input
%
% Revision 1.2  2006/05/01 19:17:24  roboos
% better support for string inputs
%
% Revision 1.1  2006/05/01 19:14:12  roboos
% first implementation as stand-alone function
%

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

