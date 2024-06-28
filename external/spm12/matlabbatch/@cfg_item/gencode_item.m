function [str, tag, cind, ccnt] = gencode_item(item, tag, tagctx, stoptag, tropts)

% function [str, tag, cind, ccnt] = gencode_item(item, tag, tagctx, stoptag, tropts)
% Generate code to recreate a generic item. This code should be suitable
% for all derived classes. Derived classes that add their own fields should
% first call this code and then add code to recreate their additional
% fields. This code does not deal with arrays of cfg_items, such a
% configuration should not exist with the current definition of a
% configuration tree.
%
% Traversal options
% struct with fields
% stopspec - match spec to stop code generation
% dflag    - (not used here)
% clvl     - current level in tree
% mlvl     - maximum level to generate - range 1 (top level only) to
%            Inf (all levels)
% cnt      - item count - used for unique tags
% mcnt     - (not evaluated here)
% Code generation stops at this item, if item matches tropts.stopspec or
% tropts.clvl > tropts.mlvl. In this case, the tag of the item is
% generated from genvarname(sprintf('%s%s', stoptag, tag), tagctx), but
% no code is generated.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode_item.m 7473 2018-11-06 10:26:44Z guillaume $

rev = '$Rev: 7473 $'; %#ok

%% Class of item
% if there are function handles in .check or .def, add their names to
% tagctx
if ~isempty(item.check) && isa(item.check, 'function_handle')
    functx = {func2str(item.check)};
else
    functx = {};
end
if ~isempty(item.def) && isa(item.def, 'function_handle')
    functx{end+1} = func2str(item.def);
end
tagctx = [tagctx functx];
% Check whether to generate code
if (tropts.clvl > tropts.mlvl || (~isempty(tropts.stopspec) && match(item, tropts.stopspec)))
    % Stopping - tag based on stoptag, and tag of item
    if isempty(tag)
        tag = genvarname(sprintf('%s%s', stoptag, item.tag), tagctx);
    else
        tag = genvarname(sprintf('%s%s', stoptag, tag), tagctx);
    end
    str = {};
    cind = [];
    ccnt = 0;
    return;
else
    if isempty(tag)
        tag = genvarname(item.tag, tagctx);
    elseif ~isletter(tag(1))
        % only modify tag if there seems to be something wrong
        % could be more specific (parse struct,cell,array subscripts)
        tag = genvarname(tag, tagctx);
    end
end
tagctx = [tagctx {tag}];
% Item count
ccnt = 1;
% Generate generic object
str{1} = '% ---------------------------------------------------------------------';
str{2} = sprintf('%% %s %s', item.tag, item.name);
str{3} = '% ---------------------------------------------------------------------';
str{4} = sprintf('%s         = %s;', tag, class(item));
% save position of class definition (needs to be overwritten if derived
% items call this generic function on their cfg_item field)
cind = 4;
%% Tag and Name
% Set basic fields
str{end+1} = sprintf('%s.tag     = ''%s'';', tag, item.tag);
% Add three spaces to name tag - this will align equal signs
% Works only because gencode does not produce subscripts for strings
str1 = gencode(item.name, sprintf('%s.name   ', tag), tagctx);
str = [str(:)' str1(:)'];
%% Val
% Generate val field
if numel(item.val) > 0 && isa(item.val{1}, 'cfg_item')
    % Traverse val{:} tree, if items are cfg_items
    cstr = {};
    % Update clvl
    ctropts = tropts;
    ctropts.clvl = ctropts.clvl + 1;
    ctropts.cnt  = ctropts.cnt + ccnt;
    ctag = cell(size(item.val));
    for k = 1:numel(item.val)
        % tags are used as variable names and need to be unique in the
        % context of this .val tag. This includes the item's tag itself
        % and the tags of its children.
        ctag{k} = genvarname(subsref(item.val{k}, substruct('.','tag')), ...
                             tagctx);
        [ccstr, ctag{k}, ccind, cccnt] = gencode_item(item.val{k}, ctag{k}, tagctx, ...
                                              stoptag, ctropts);
        if ~isempty(ccstr)
            % Child has returned code
            cstr = [cstr(:)' ccstr(:)'];
            ccnt = ccnt + cccnt;
            ctropts.cnt = ctropts.cnt + cccnt;
            tagctx = [tagctx ctag(k)];
        end
    end
    % Update position of class definition
    cind = cind+numel(cstr);
    % Prepend code of children
    str = [cstr(:)' str(:)'];
    str{end+1} = sprintf('%s.val     = {%s};' ,tag, sprintf('%s ', ctag{:}));
elseif numel(item.val) > 0 && ~isa(item.val{1}, 'cfg_item')
    % Check .def field. Generate code for .val only, if no defaults
    % defined or value is different from defaults.
    if exist('isequalwithequalnans','builtin')
        iseqn = isequalwithequalnans(feval(item.def), item.val{1});
    else
        iseqn = isequaln(feval(item.def), item.val{1});
    end
    if isempty(item.def) || ~iseqn
        str1 = gencode(item.val, sprintf('%s.val', tag), tagctx);
        str = [str(:)' str1(:)'];
    end
end
%% Check
% Generate check field
if ~isempty(item.check)
    % Add two spaces to check tag - this will align equal signs
    % Works only because gencode does not produce subscripts for function
    % strings
    str1 = gencode(item.check, sprintf('%s.check  ', tag), tagctx);
    str = [str(:)' str1(:)'];
end    
%% Rewrite job
% Generate rewrite_job field
if ~isempty(item.rewrite_job)
    str1 = gencode(item.rewrite_job, sprintf('%s.rewrite_job', tag), tagctx);
    str = [str(:)' str1(:)'];
end    
%% Help
% Generate help field
if numel(item.help) > 0
    % Add three spaces to help tag - this will align equal signs
    % Works only because gencode does not produce subscripts for cellstrings
    str1 = gencode(item.help, sprintf('%s.help   ', tag), tagctx);
    str = [str(:)' str1(:)'];
end
%% Def
% Generate def field
if ~isempty(item.def)
    % Add four spaces to def tag - this will align equal signs
    % Works only because gencode does not produce subscripts for function
    % strings
    str1 = gencode(item.def, sprintf('%s.def    ', tag), tagctx);
    str = [str(:)' str1(:)'];
end
