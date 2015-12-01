function tree = move(tree,uida, uidb)
% XMLTREE/MOVE Move (move a subtree inside a tree from A to B)
% 
% tree   - XMLTree object
% uida   - initial position of the subtree
% uidb   - parent of the final position of the subtree
%__________________________________________________________________________
%
% Move a subtree inside a tree from A to B.
% The tree parameter must be in input AND in output.
%__________________________________________________________________________
% Copyright (C) 2002-2015  http://www.artefact.tk/

% Guillaume Flandin
% $Id$


%error(nargchk(3,3,nargin));

p = tree.tree{uida}.parent;
tree.tree{p}.contents(find(tree.tree{p}.contents==uida)) = [];
tree.tree(uidb).contents = [tree.tree(uidb).contents uida];
