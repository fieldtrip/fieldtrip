function item = setval(item, val, dflag)

% function item = setval(item, val, dflag)
% set item.val{1} to val. Validity checks are performed through subsasgn.
% If val == {}, set item.val to {}. 
% If dflag is true, and item.def is not empty, set the default setting for
% this item instead by calling feval(item.def{:}, val). If val == {}, use
% the string '<UNDEFINED>' as in a harvested tree. If dflag is true, but
% no item.def defined, set item.val{1} instead.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: setval.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

if iscell(val) && isempty(val)
    if dflag
        if ~isempty(item.def)
            try
                feval(item.def, {'<UNDEFINED>'});
            catch
                cfg_message('matlabbatch:setval:defaults', ...
                            '%s: unable to set default value.', ...
                            subsasgn_checkstr(item, substruct('.','val')));
            end
        else
            item = subsasgn(item, substruct('.','val'), {});
        end
    else
        item = subsasgn(item, substruct('.','val'), {});
    end
else
    if dflag
        [sts, val1] = subsasgn_check(item, substruct('.','val'), {val});
        if sts
            if ~isempty(item.def)
                try
                    feval(item.def, val1);
                catch
                    cfg_message('matlabbatch:setval:defaults', ...
                                '%s: unable to set default value.', ...
                                subsasgn_checkstr(item, substruct('.','val')));
                end
            else
                item = subsasgn(item, substruct('.','val'), {val});
            end
        end
    else
        item = subsasgn(item, substruct('.','val'), {val});
    end
end
