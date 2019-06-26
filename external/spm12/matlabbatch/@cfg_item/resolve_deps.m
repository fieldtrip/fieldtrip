function [val, sts] = resolve_deps(item, cj)

% function [val, sts] = resolve_deps(item, cj)
% Resolve dependencies for an cfg item. This is a generic function that
% returns the contents of item.val{1} if it is an array of cfg_deps. If
% there is more than one dependency, they will be resolved in order of
% appearance. The returned val will be the concatenation of the values of
% all dependencies. A warning will be issued if this concatenation fails
% (which would happen if resolved dependencies contain incompatible
% values).
% If any of the dependencies cannot be resolved, val will be empty and sts
% false.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: resolve_deps.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

val1 = cell(size(item.val{1}));
for k = 1:numel(item.val{1})
    % Outputs are stored in .jout field of cfg_exbranch, which is
    % not included in .src_exbranch substruct
    out = subsref(cj, [item.val{1}(k).src_exbranch, ...
                            substruct('.','jout')]);
    sts = ~isa(out,'cfg_inv_out');
    if ~sts
        % dependency not yet computed, fail silently
        val = [];
        return;
    end
    try
        val1{k} = subsref(out, item.val{1}(k).src_output);
    catch 
        % dependency can't be resolved, even though it should be there
        l = lasterror; %#ok
        % display source output to diagnose problems
        val1{k} = out;
        dstr = disp_deps(item, val1);
        cfg_message('matlabbatch:resolve_deps:missing', ...
                'Dependency source available, but is missing required output.\n%s',...
                l.message);
        cfg_message('matlabbatch:resolve_deps:missing', '%s\n', dstr{:});
        val = [];
        sts = false;
        return;
    end
end
if sts
    % All items resolved, try concatenation
    try
        % try concatenation along 1st dim
        val = cat(1, val1{:});
    catch 
        % try concatenation along 2nd dim
        try
            val = cat(2, val1{:});
        catch 
            % all concatenations failed, display warning
            l = lasterror; %#ok
            dstr = disp_deps(item, val1);
            cfg_message('matlabbatch:resolve_deps:incompatible',...
                        'Dependencies resolved, but incompatible values.\n%s', ...
                        l.message);
            cfg_message('matlabbatch:resolve_deps:incompatible', '%s\n', dstr{:});
            % reset val and sts
            val = [];
            sts = false;
            return;
        end
    end
end
% all collected, check subsasgn validity
if sts
    % subsasgn_check only accepts single subscripts
    [sts, val] = subsasgn_check(item, substruct('.','val'), {val});
end;
if sts
    % dereference val after subsasgn_check
    val = val{1};
else
    dstr = disp_deps(item, val1);
    dstr = [{'Dependencies resolved, but not suitable for this item.'}, ...
        dstr(:)'];
    cfg_message('matlabbatch:subsasgn:val',...
        '%s\n', dstr{:});
    return;
end

function dstr = disp_deps(item, val1) %#ok<INUSD>
dstr = cell(numel(item.val{1})+1,1);
dstr{1} = sprintf('In item %s:', subsref(item, substruct('.','name')));
for k = 1:numel(item.val{1})
    substr = gencode_substruct(item.val{1}(k).src_output);
    dstr{k+1} = sprintf('Dependency %d: %s (out%s)\n%s', ...
            k, item.val{1}(k).sname, substr{1}, evalc('disp(val1{k})'));
end
