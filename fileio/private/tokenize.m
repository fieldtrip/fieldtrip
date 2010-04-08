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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

