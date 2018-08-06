function s = char(tree)
% XMLTREE/CHAR Converter function from XMLTree to a description string
% FORMAT s = char(tree)
%
% tree - XMLTree object
% s    - a description string of an XMLTree
%__________________________________________________________________________
%
% Return a string describing the XMLTree:
%               'XMLTree object (x nodes) [filename]'
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: char.m 4460 2011-09-05 14:52:16Z guillaume $


s = strcat('XMLTree object (',num2str(length(tree)),' nodes) [',getfilename(tree),']');
