function varargout = add(tree,uid,type,parameter)
% XMLTREE/ADD Method (add childs to elements of an XML Tree)
% FORMAT vararout = add(tree,uid,type,parameter)
% 
% tree      - XMLTree object
% uid       - array of uid's
% type      - 'element', 'chardata', 'cdata', 'pi' or 'comment'
% parameter - property name (a character array unless type='pi' for
%             which parameter=struct('target','','value',''))
%
% new_uid   - UID's of the newly created nodes
%
%        tree = add(tree,uid,type,parameter);
%        [tree, new_uid] = add(tree,uid,type,parameter);
%__________________________________________________________________________
%
% Add a node (element, chardata, cdata, pi or comment) in the XML Tree.
% It adds a child to the element whose UID is iud.
% Use attributes({'set','get','add','del','length'},...) function to 
% deal with the attributes of an element node (initialized empty).
% The tree parameter must be in input AND in output.
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: add.m 4460 2011-09-05 14:52:16Z guillaume $


if ~isa(uid,'double')
    error('[XMLTree] UID must be a double array.');
end
if ~ischar(type)
    error('[XMLTree] TYPE must be a valid item type.');
end
if strcmp(type,'pi')
    if ~isfield(parameter,'target') || ~isfield(parameter,'value') || ...
       ~ischar(parameter.target) || ~ischar(parameter.value)
        error(['[XMLTree] For a Processing Instruction, ',...
                        'PARAMETER must be a struct.']);
    end
elseif ~ischar(parameter)
    error('[XMLTree] PARAMETER must be a string.');
end

if nargout == 2
    l = length(tree.tree);
    varargout{2} = (l+1):(l+numel(uid));
end

for i=1:numel(uid)
    if uid(i)<1 || uid(i)>length(tree.tree)
        error('[XMLTree] Invalid UID.');
    end
    if ~strcmp(tree.tree{uid(i)}.type,'element')
        error('[XMLTree] Cannot add a child to a non-element node.');
    end
    l = length(tree.tree);
    switch type
        case 'element'
            tree.tree{l+1} = struct('type','element',...
                                    'name',parameter,...
                                    'attributes',[],...
                                    'contents',[],...
                                    'parent',[],...
                                    'uid',l+1);
        case 'chardata'
            tree.tree{l+1} = struct('type','chardata',...
                                    'value',parameter,...
                                    'parent',[],...
                                    'uid',l+1);
        case  'cdata'
            tree.tree{l+1} = struct('type','cdata',...
                                    'value',parameter,...
                                    'parent',[],...
                                    'uid',l+1);
        case 'pi'
            tree.tree{l+1} = struct('type','pi',...
                                    'target',parameter.target,...
                                    'value',parameter.value,...
                                    'parent',[],...
                                    'uid',l+1);
        case 'comment'
            tree.tree{l+1} = struct('type','comment',...
                                    'value',parameter,...
                                    'parent',[],...
                                    'uid',l+1);
        otherwise
            error(sprintf('[XMLTree] %s: unknown item type.',type));
    end
    tree.tree{uid(i)}.contents = [tree.tree{uid(i)}.contents l+1];
    tree.tree{l+1}.parent = uid(i);
end

varargout{1} = tree;
