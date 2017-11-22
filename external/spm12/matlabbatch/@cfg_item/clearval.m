function item = clearval(item, dflag)

% function item = clearval(item, dflag)
% This is a generic function to clear the contents of the val field of a
% cfg_item. It is usable for all leaf cfg_item classes (cfg_entry,
% cfg_files, cfg_menu).
% dflag is ignored for leaf entries - the val field is cleared
% unconditionally.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: clearval.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

item.val = {};
