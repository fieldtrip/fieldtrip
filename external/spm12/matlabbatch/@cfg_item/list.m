function [id, stop, val] = list(item, spec, tropts, fn)

% function [id, stop, val] = list(item, spec, tropts, fn)
% This function searches the cfg tree for certain entries.
%
% [id stop val] = list(item, spec, tropts[, fieldname])
% Find items in a cfg tree rooted at item that match a specification spec.
% By default, the filled configuration tree is searched (i.e. the
% val-branches of cfg_repeat and cfg_choice nodes). 
% See MATCH for help about spec data structure.
%
% Traversal options
% struct with fields
% stopspec - match spec to stop traversal
% dflag    - traverse val or values tree
% clvl     - current level in tree
% mlvl     - maximum level to traverse - range 1 (top level only) to
%            Inf (all levels)
% cnt      - #items found so far
% mcnt    - max #items to find
% List will stop descending into subtrees if one of the conditions
% following conditions are met: item matches stopspec, clvl >= mlvl, cnt >=
% mcnt. Flag stop is true for nodes where traversal has stopped
% (i.e. items where tropts has stopped further traversal).
%
% A cell list of subsref ids to matching nodes will be returned. The id of
% this node is returned before the id of its matching children.
% If the root node of the tree matches, the first id returned will be an
% empty substruct.
% If a cell list of fieldnames is given, then the contents of these fields
% will be returned in the cell array val. If one of the fields does not
% exist, a cell with an empty entry will be returned.
% There are five pseudo-fieldnames which allow to obtain information useful
% to build e.g. a user interface for cfg trees:
% 'class' - returns the class of the current item
% 'level' - returns the level in the tree. Since data is collected
%           pre-order, children are listed after their parents. Identical
%           levels of subsequent nodes denote siblings, whereas decreasing
%           levels of subsequent nodes denote siblings of the parent node.
% 'all_set' - return all_set status of subtree rooted at item, regardless
%             whether list will descend into it or not
% 'all_set_item' - return all_set_item status of current node (i.e. whether
%                  all integrity conditions for this node are fulfilled)
%                  For in-tree nodes this can be different from all_set.
% 'showdoc' - calls showdoc to display the help text and option hints for
%             the current item.
% This code is the generic list function, suitable for all cfg_leaf items.
% To ensure that the correct val (val{1}, dependency or default value)
% is listed, the val field is treated in a special way.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: list.m 3469 2009-10-16 08:43:15Z volkmar $

rev = '$Rev: 3469 $'; %#ok

if match(item, spec)
    id = {struct('type', {}, 'subs', {})};
    stop = false;
    if nargin > 3
        val = cell(1,numel(fn));
        for k = 1:numel(fn)
            switch fn{k}
                case 'class'
                    val{k} = {class(item)};
                case 'level'
                    val{k} = {tropts.clvl};
                case 'all_set'
                    val{k} = {all_set(item)};
                case 'all_set_item'
                    val{k} = {all_set_item(item)};
                case 'showdoc'
                    val{k} = {showdoc(item,'')};
                case 'val'
                    % special case for val item
                    if tropts.dflag
                        if isempty(item.def)
                            dval = item.val;
                            if isa(dval, 'cfg_dep')
                                dval = [];
                            end;
                            val{k} = {dval};
                        else
                            val{k} = {{item.def({})}};
                        end
                    else
                         val{k} = {subsref(item, substruct('.', fn{k}))};
                    end
                case fieldnames(item)
                    val{k} = {subsref(item, substruct('.', fn{k}))};
                otherwise
                    val{k} = {{}};
            end;
        end;
    else
        val = {};
    end;
else
    id = {};
    stop = [];
    if nargin > 3
        val = cell(size(fn));
        [val{:}] = deal({});        
    else
        val = {};
    end;
end;
