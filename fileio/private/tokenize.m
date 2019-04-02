function [tok] = tokenize(str, sep, rep)

% TOKENIZE cuts a string into pieces, returning the pieces in a cell-array
%
% Use as
%   t = tokenize(str)
%   t = tokenize(str, sep)
%   t = tokenize(str, sep, rep)
% where
%   str = the string that you want to cut into pieces
%   sep = the separator at which to cut (default is whitespace)
%   rep = whether to treat repeating separator characters as one (default is false)
%
% With the optional boolean flag "rep" you can specify whether repeated
% separator characters should be squeezed together (e.g. multiple
% spaces between two words). The default is rep=1, i.e. repeated
% separators are treated as one.
%
% See also STRTOK, TEXTSCAN

% Copyright (C) 2003-2010, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

str = str(:)';

if nargin<2
  sep = [9:13 32]; % White space characters
end

if nargin<3
  rep = false;
end

current_argin = {str, sep, rep};
if isequal(current_argin, previous_argin)
  % don't do the processing again, but return the previous values from cache
  tok = previous_argout;
  return
end

if numel(sep)==1
  f = find(str==sep);
else
  f = find(ismember(str, sep));
end
f = [0, f, length(str)+1];

tok = cell(1, length(f)-1);
for i=1:(length(f)-1)
  tok{i} = str((f(i)+1):(f(i+1)-1));
end

if rep
  % remove empty cells, which occur if the separator is repeated (e.g. multiple spaces)
  tok(cellfun('isempty', tok))=[];
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout  = tok;
previous_argin  = current_argin;
previous_argout = current_argout;

