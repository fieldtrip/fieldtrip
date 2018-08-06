function child = children(tree,uid)
% XMLTREE/CHILDREN Return children's UIDs of node uid
% FORMAT child = children(tree,uid)
%
% tree   - a tree
% uid    - uid of the element
% child  - array of the UIDs of children of node uid
%__________________________________________________________________________
%
% Return UID's of children of node uid
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: children.m 4460 2011-09-05 14:52:16Z guillaume $

%error(nargchk(2,2,nargin));

child = [];
uid = uid(:);
l = length(tree);
for i=1:length(uid)
    if uid(i) > 0 && uid(i) <= l
        if strcmp(tree.tree{uid(i)}.type,'element')
            child = [child tree.tree{uid(i)}.contents];
        end
    else
        error('[XMLTree] Invalid UID.');
    end
end
if isempty(child), child = []; end
