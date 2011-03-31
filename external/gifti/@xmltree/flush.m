function tree = flush(tree,uid)
% XMLTREE/FLUSH Flush (Clear a subtree given its UID)
% 
% tree      - XMLTree object
% uid       - array of UID's of subtrees to be cleared
%             Default is root
%_______________________________________________________________________
%
% Clear a subtree given its UID (remove all the leaves of the tree)
% The tree parameter must be in input AND in output
%_______________________________________________________________________
% Copyright (C) 2002-2008  http://www.artefact.tk/

% Guillaume Flandin <guillaume@artefact.tk>
% $Id$

error(nargchk(1,2,nargin));

if nargin == 1,
    uid = root(tree);
end

uid = uid(:);
for i=1:length(uid)
     tree = sub_flush(tree,uid(i));
end

%=======================================================================
function tree = sub_flush(tree,uid)
    if isfield(tree.tree{uid},'contents')
        % contents is parsed in reverse order because each child is
        % deleted and the contents vector is then eventually reduced
        for i=length(tree.tree{uid}.contents):-1:1
            tree = sub_flush(tree,tree.tree{uid}.contents(i));
        end
    end
    if strcmp(tree.tree{uid}.type,'chardata') ||...
        strcmp(tree.tree{uid}.type,'pi') ||...
        strcmp(tree.tree{uid}.type,'cdata') ||...
        strcmp(tree.tree{uid}.type,'comment')
        tree = delete(tree,uid);
    end
