function varargout = attributes(varargin)
% XMLTREE/ATTRIBUTES Method (handle attributes of an element node)
% FORMAT varargout = attributes(varargin)
% 
% tree    - XMLTree object
% method  - 'set','get','add','del' or 'length'
% uid     - the UID of an element node
% n       - indice of the attribute
% key     - string key="..."
% val     - string ...="val"
% attr    - cell array of struct(key,val) or just struct(key,val)
% l       - number of attributes of the element node uid
%
%     tree = attributes(tree,'set',uid,n,key,val)
%     attr = attributes(tree,'get',uid[,n])
%     tree = attributes(tree,'add',uid,key,val)
%     tree = attributes(tree,'del',uid[,n])
%     l    = attributes(tree,'length',uid)
%__________________________________________________________________________
%
% Handle attributes of an element node.
% The tree parameter must be in input AND in output for 'set', 'add' and
% 'del' methods.
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: attributes.m 4460 2011-09-05 14:52:16Z guillaume $


%error(nargchk(3,6,nargin));

tree = varargin{1};
if ~ischar(varargin{2}) || ...
   ~any(strcmp(varargin{2},{'set','get','add','del','length'}))
    error('[XMLTree] Unknown method.');
end
uid = varargin{3};
if ~isa(uid,'double') || any(uid>length(tree)) || any(uid<1)
    error('[XMLTree] UID must be a positive integer scalar.');
end

if ~strcmp(tree.tree{uid}.type,'element')
    error('[XMLTree] This node has no attributes.');
end

switch varargin{2}
    case 'set'
        %error(nargchk(6,6,nargin));
        if ~isa(varargin{4},'double') || ...
           any(varargin{4}>length(tree.tree{uid}.attributes)) || ...
           any(varargin{4}<1)
            error('[XMLTree] Invalid attribute indice.');
        end
        ind = varargin{4};
        tree.tree{uid}.attributes{ind} = struct('key',varargin{5},'val',varargin{6});
        varargout{1} = tree;
    case 'get'
        %error(nargchk(3,4,nargin));
        if nargin == 4
            if ischar(varargin{4})
                for i=1:length(tree.tree{uid}.attributes)
                    if strcmp(varargin{4},tree.tree{uid}.attributes{i}.key)
                        varargout{1} = tree.tree{uid}.attributes{i}.val;
                        return;
                    end
                end
                varargout{1} = [];
            elseif ~isa(varargin{4},'double') || ...
               any(varargin{4}>length(tree.tree{uid}.attributes)) || ...
               any(varargin{4}<1)
                error('[XMLTree] Invalid attribute indice.');
            else
                if length(varargin{4}) == 1
                    varargout{1} = tree.tree{uid}.attributes{varargin{4}(1)};
                else
                    varargout{1} = {};
                    for i=1:length(varargin{4})
                        varargout{1}{i} = tree.tree{uid}.attributes{varargin{4}(i)};
                    end
                end
            end
        else
            if length(tree.tree{uid}.attributes) == 1
                varargout{1} = tree.tree{uid}.attributes{1};
            else
                varargout{1} = {};
                for i=1:length(tree.tree{uid}.attributes)
                    varargout{1}{i} = tree.tree{uid}.attributes{i};
                end
            end
        end
    case 'add'
        %error(nargchk(5,5,nargin));
        ind = length(tree.tree{uid}.attributes) + 1;
        tree.tree{uid}.attributes{ind} = struct('key',varargin{4},'val',varargin{5});
        varargout{1} = tree;
    case 'del'
        %error(nargchk(3,4,nargin));
        if nargin == 4
            if ~isa(varargin{4},'double') || ...
              any(varargin{4}>length(tree.tree{uid}.attributes)) || ...
              any(varargin{4}<1)
                error('[XMLTree] Invalid attribute indice.');
            end
            ind = varargin{4};
            tree.tree{uid}.attributes(ind) = [];
        else
            tree.tree{uid}.attributes = [];
        end
        varargout{1} = tree;
    case 'length'
        %error(nargchk(3,3,nargin));
        varargout{1} = length(tree.tree{uid}.attributes);
    otherwise
        error('[XMLTree] Unknown method.');
end
