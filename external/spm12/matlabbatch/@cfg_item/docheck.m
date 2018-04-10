function chk = docheck(item, val)

% function chk = docheck(item, val)
% Run item specific check function, if present.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: docheck.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

chk = true;

if ~isempty(item.check) && all_set(item)
    % Call content dependent check function, if it is present and
    % dependencies are to be resolved
    try
        cstr = feval(item.check, val);
        chk  = isempty(cstr);
    catch
        cstr = sprintf('Check function ''%s'' failed.', ...
                       func2str(item.check));
        chk  = false;
    end;
    if ~chk
        if iscellstr(cstr)
            cstr = sprintf('%s\n', cstr{:});
        end;
        cfg_message('matlabbatch:harvest:check', ...
                'Contents of ''%s'' does not meet check criteria:\n''%s''', ...
                gettag(item), cstr);
    end;
end;