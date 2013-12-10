function display(tree)
% XMLTREE/DISPLAY Command window display of an XMLTree
% FORMAT display(tree)
% 
% tree - XMLTree object
%__________________________________________________________________________
%
% This method is called when the semicolon is not used to terminate a
% statement which returns an XMLTree.
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id$

disp(' ');
disp([inputname(1),' = ']);
disp(' ');
for i=1:numel(tree)
    disp([blanks(length(inputname(1))+3) char(tree(i))]);
end
disp(' ');
