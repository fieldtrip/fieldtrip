function tree = set(tree,uid, parameter, value)
% XMLTREE/SET Method (set object properties)
% FORMAT tree = set(tree,uid,parameter,value)
% 
% tree      - XMLTree object
% uid       - array (or cell) of uid's
% parameter - property name
% value     - property value
%__________________________________________________________________________
%
% Set object properties given its uid and pairs parameter/value
% The tree parameter must be in input AND in output
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: set.m 4460 2011-09-05 14:52:16Z guillaume $


%error(nargchk(4,4,nargin));

if iscell(uid), uid = [uid{:}]; else uid = uid(:); end

for i=1:length(uid)
    tree.tree{uid(i)} = builtin('subsasgn', tree.tree{uid(i)}, struct('type','.','subs',parameter), value);
    %tree.tree{uid(i)} = setfield(tree.tree{uid(i)},parameter,value);
end
