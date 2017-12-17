function [item, inputs] = fillvals(item, inputs, infcn)

% function [item, inputs] = fillvals(item, inputs, infcn)
% If ~all_set_item, try to set item.val{1} to inputs{1}. Validity checks
% are performed through subsasgn. If inputs{1} is not suitable for this
% item, it is discarded. If infcn is a function handle,
% [val sts] = infcn(item) 
% will be called to obtain a value for this item. This call will be
% repeated until either val can be assigned to item or sts is true. val
% should be a cell array with 1 item and val{1} the value that is to be
% assigned to item.val{1}.
% This function is suitable for all cfg_leaf items.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: fillvals.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

if ~all_set_item(item)
    if ~isempty(inputs)
        item = setval(item, inputs{1}, false);
        inputs = inputs(2:end);
    end;
    if ~all_set_item(item) && ~isempty(infcn) && subsasgn_check_funhandle(infcn)
        sts = false;
        while ~sts && ~all_set_item(item)
            [val, sts] = feval(infcn, item);
            if sts
                item = setval(item, val, false);
            end;
        end;
    end;
end;
