function [tok] = tokenize(str, sep, rep)

% TOKENIZE cuts a string into pieces, returning the pieces in a cell array
%
% Use as
%   t = tokenize(str)
%   t = tokenize(str, sep)
%   t = tokenize(str, sep, rep)
% where
%   str    = the string that you want to cut into pieces
%   sep    = the separator at which to cut (default is whitespace)
%   rep    = whether to treat repeating seperator characters as one (default is false)
%
% With the optional boolean flag "rep" you can specify whether repeated
% seperator characters should be squeezed together (e.g. multiple
% spaces between two words). The default is rep=1, i.e. repeated
% seperators are treated as one.
%
% See also STRTOK, TEXTSCAN

% Copyright (C) 2003-2008, Robert Oostenveld
%
% $Log: tokenize.m,v $
% Revision 1.3  2008/11/20 15:46:04  roboos
% fixed bug in the default for the 3rd argument (whether repeated seperators should be treated as one). This bug caused some trouble with reading brainvision header files over the last two weeks.
%
% Revision 1.2  2008/11/14 07:37:05  roboos
% use whitespace if no seperator is specified
%
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.5  2006/05/10 07:15:22  roboos
% allow for squeezing multiple separators into one
%
% Revision 1.4  2006/01/30 14:56:41  roboos
% fixed another stupid bug in previous cvs commit
%
% Revision 1.3  2006/01/30 14:55:44  roboos
% fixed stupid bug in previous cvs commit
%
% Revision 1.2  2006/01/30 13:38:18  roboos
% replaced dependency on strtok by more simple code
% changed from dos to unix
%
% Revision 1.1  2005/05/23 13:47:51  roboos
% old implementation, new addition to CVS for fieldtrip release
%

if nargin<2
  sep = [9:13 32]; % White space characters
end

if nargin<3
  rep = false;
end

tok = {};
f = find(ismember(str, sep));
f = [0, f, length(str)+1];
for i=1:(length(f)-1)
  tok{i} = str((f(i)+1):(f(i+1)-1));
end

if rep
  % remove empty cells, which occur if the separator is repeated (e.g. multiple spaces)
  tok(cellfun('isempty', tok))=[];
end

