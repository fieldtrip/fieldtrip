function item = subsasgn_job(item, subs, val)

% function item = subsasgn_job(item, subs, val)
% Treat a subscript reference as a reference in a job structure instead
% of a cfg_item structure. This generic cfg_item method treats subs as a
% subscript reference into item.val{1}. This is suitable for all cfg_leaf
% items.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_job.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok

item = subsasgn(item, [substruct('.','val','{}',{1}) subs], val);
