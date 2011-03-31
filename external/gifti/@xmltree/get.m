function value = get(tree,uid,parameter)
% XMLTREE/GET Get Method (get object properties)
% FORMAT value = get(tree,uid,parameter)
% 
% tree      - XMLTree object
% uid       - array of uid's
% parameter - property name
% value     - property value
%_______________________________________________________________________
%
% Get object properties of a tree given their UIDs.
%_______________________________________________________________________
% Copyright (C) 2002-2008  http://www.artefact.tk/

% Guillaume Flandin <guillaume@artefact.tk>
% $Id$

error(nargchk(2,3,nargin));

value = cell(size(uid));
uid = uid(:);
if nargin==2
    for i=1:length(uid)
        if uid(i)<1 || uid(i)>length(tree.tree)
            error('[XMLTree] Invalid UID.');
        end
        % According to the type of the node, return only some parameters
        % Need changes...
        value{i} = tree.tree{uid(i)};
    end
else
    for i=1:length(uid)
        try
            value{i} = subsref(tree.tree{uid(i)}, struct('type','.','subs',parameter));
        catch
            error(sprintf('[XMLTree] Parameter %s not found.',parameter));
        end
    end 
end
if length(value)==1
    value = value{1};
end  
