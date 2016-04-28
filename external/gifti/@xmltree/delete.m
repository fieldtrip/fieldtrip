function tree = delete(tree,uid)
% XMLTREE/DELETE Delete (delete a subtree given its UID)
% 
% tree      - XMLTree object
% uid       - array of UID's of subtrees to be deleted
%__________________________________________________________________________
%
% Delete a subtree given its UID
% The tree parameter must be in input AND in output
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: delete.m 4460 2011-09-05 14:52:16Z guillaume $

%error(nargchk(2,2,nargin));

uid = uid(:);
for i=1:length(uid)
    if uid(i)==1
        warning('[XMLTree] Cannot delete root element.');
    else
        p = tree.tree{uid(i)}.parent;
        tree = sub_delete(tree,uid(i));
        tree.tree{p}.contents(find(tree.tree{p}.contents==uid(i))) = [];
    end
end

%==========================================================================
function tree = sub_delete(tree,uid)
    if isfield(tree.tree{uid},'contents')
        for i=1:length(tree.tree{uid}.contents)
            tree = sub_delete(tree,tree.tree{uid}.contents(i));
        end
    end
    tree.tree{uid} = struct('type','deleted');
