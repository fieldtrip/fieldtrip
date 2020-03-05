function ok = all_leafs(item)

% function ok = all_leafs(item)
% Generic all_leafs function that returns true. This is suitable for all
% leaf items. No content specific checks are performed.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: all_leafs.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok
% do not check anything else than item class
ok = isa(item, 'cfg_leaf');

