function ok = all_set(item)

% function ok = all_set(item)
% Generic all_set function - checks whether item.val is not empty. No
% checks based on the content of item.val are performed here.
% Content checking is done in the following places:
% * context-insensitive checks based on configuration specifications
%   are performed during subsasgn/setval. This will happen during user
%   input or while resolving dependencies during harvest. 
% * context sensitive checks by a configuration .check function are
%   performed during harvest after all dependencies are resolved.
% This function is suitable for all leaf configuration items.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: all_set.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok
% by default, do not check input size/type etc. this is already done in
% subsasgn
ok = ~isempty(item.val)||item.hidden;
