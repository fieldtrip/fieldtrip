function tree = copy(tree,subuid,uid)
% XMLTREE/COPY Copy Method (copy a subtree in another branch)
% FORMAT tree = copy(tree,subuid,uid)
% 
% tree      - XMLTree object
% subuid    - UID of the subtree to copy
% uid       - UID of the element where the subtree must be duplicated
%__________________________________________________________________________
%
% Copy a subtree to another branch.
% The tree parameter must be in input AND in output.
%__________________________________________________________________________
% Copyright (C) 2002-2015  http://www.artefact.tk/

% Guillaume Flandin
% $Id: copy.m 6480 2015-06-13 01:08:30Z guillaume $


%error(nargchk(2,3,nargin));

if nargin == 2
    uid = parent(tree,subuid);
end

l = length(tree);
tree = sub_copy(tree,subuid,uid);
tree.tree{uid}.contents = [tree.tree{uid}.contents l+1];

% to have the copy next to the original and not at the end?
%  contents = get(tree,parent,'contents');
%  i = find(contents==uid);
%  tree = set(tree,parent,'contents',[contents(1:i) l+1 contents(i+1:end)]);

%==========================================================================
function tree = sub_copy(tree,uid,p)

    l = length(tree);
    tree.tree{l+1} = tree.tree{uid};
    tree.tree{l+1}.uid = l+1;
    tree.tree{l+1}.parent = p;
    tree.tree{l+1}.contents = [];
    if isfield(tree.tree{uid},'contents')
        contents = get(tree,uid,'contents');
        m = length(tree);
        for i=1:length(contents)
            tree.tree{l+1}.contents = [tree.tree{l+1}.contents m+1];
            tree = sub_copy(tree,contents(i),l+1);
            m = length(tree);
        end
    end
