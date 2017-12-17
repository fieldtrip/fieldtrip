function [item, sts] = expand(item, eflag, tropts)

% function [item, sts] = expand(item, eflag, tropts)
% Set/query expanded flag of item depending on eflag:
% -1 - do not force eflag to any state, only child state will be inherited
%  0 - collapse
%  1 - expand val unconditionally
%  2 - expand metadata unconditionally
%  3 - expand val, if it is not set
% Return status is (expanded > 0), i.e. if expanded, then no additional
% info about expansion level or expansion reason is returned and parent
% nodes are set to expanded = 1.
%
% Traversal options
% struct with fields
% stopspec - match spec to stop traversal
% dflag    - traverse val or values tree
% clvl     - current level in tree
% mlvl     - maximum level to traverse - range 1 (top level only) to
%            Inf (all levels)
% cnt (not set here)
% mcnt (not evaluated here)
% Traversal options are used here to control which items should be forced
% to expand/unexpand. Traversal continues to child items, even if level or
% stopspec criteria are met, but with an eflag of -1 (i.e. only 'expanded'
% status is queried, but not changed).
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: expand.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

% Set expanded based on eflag in item
if eflag >= 0 && eflag <= 2
    item.expanded = eflag;
end;
if eflag == 3
    if ~all_set(item)
        item.expanded = eflag;
    else
        item.expanded = 0;
    end;
end;
sts = item.expanded > 0;
