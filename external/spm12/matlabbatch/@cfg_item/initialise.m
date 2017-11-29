function item = initialise(item, val, dflag)

% function item = initialise(item, val, dflag)
% This is a generic initialisation function to insert values into the val
% field of a cfg_item. It is usable for all leaf cfg_item classes
% (cfg_entry, cfg_files, cfg_menu). Assignment checks are done through
% subsasgn.
% dflag is ignored for leaf entries - the passed value is assigned
% unconditionally.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: initialise.m 4898 2012-09-05 13:40:16Z volkmar $

rev = '$Rev: 4898 $'; %#ok

d1flag = false;
if ischar(val) && strcmp(val, '<DEFAULTS>')
    if isempty(item.def)
        % don't change anything if no .def field
        return;
    else
        try
            val = feval(item.def, {});
            d1flag = true;
        catch
            cfg_message('matlabbatch:initialise:defaults', ...
                    '%s: No matching defaults value found.', ...
                    subsasgn_checkstr(item,substruct('.','val')));
            val = '<UNDEFINED>';
        end;
    end;
end;

subs = substruct('.', 'val', '{}', {1});
if ischar(val) && strcmp(val, '<UNDEFINED>') % val may be <UNDEFINED>
    item.val = {};
else
    % Modify job before initialisation
    if ~d1flag && ~dflag && ~isempty(item.rewrite_job)
        val = feval(item.rewrite_job, val);
    end
    item = subsasgn(item, subs, val);
end;