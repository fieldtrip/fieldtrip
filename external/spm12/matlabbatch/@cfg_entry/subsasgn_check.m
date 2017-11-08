function [sts, val] = subsasgn_check(item,subs,val)

% function [sts, val] = subsasgn_check(item,subs,val)
% Perform validity checks for cfg_entry inputs. Does not yet support
% evaluation of inputs.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 5946 2014-04-10 13:09:34Z volkmar $

rev = '$Rev: 5946 $'; %#ok

sts = true;
switch subs(1).subs
    case {'num'}
        % special num treatment - num does describe the dimensions of
        % input in cfg_entry items, not a min/max number
        sts = isnumeric(val) && (isempty(val) || numel(val)>=2 && all(val(:) >= 0));
        if ~sts
             cfg_message('matlabbatch:check:num', ...
                     '%s: Value must be empty or a vector of non-negative numbers with at least 2 elements', ...
                     subsasgn_checkstr(item,subs));
        end
    case {'val'}
        % perform validity checks - subsasgn_check should be called with
        % a cell containing one item
        if ~iscell(val)
            cfg_message('matlabbatch:checkval', ...
                    '%s: Value must be a cell.', subsasgn_checkstr(item,subs));
            sts = false;
            return;
        end
        if isempty(val)
            val = {};
        else
            % check whether val{1} is a valid element
            [sts, vtmp] = valcheck(item,val{1});
            val = {vtmp};
        end
    case {'strtype'}
        strtypes = {'s','e','f','n','w','i','r','c','x','p'};
        sts = isempty(val) || (ischar(val) && ...
                               any(strcmp(val, strtypes)));
        if ~sts
            cfg_message('matlabbatch:check:strtype', ...
                    '%s: Value must be a valid strtype.', subsasgn_checkstr(item,subs));
        end
end

function [sts, val] = valcheck(item,val)
% taken from spm_jobman/stringval
% spm_eeval goes into GUI
sts = true;
% check for reserved words
if ischar(val) && size(val, 1) == 1 && any(strcmp(val, {'<UNDEFINED>','<DEFAULTS>'}))
    cfg_message('matlabbatch:checkval', ...
            ['%s: Item must not be one of the reserved words ''<UNDEFINED>'' ' ...
             'or ''<DEFAULTS>''.'], subsasgn_checkstr(item,substruct('.','val')));
    sts = false;
    return;
end
if isa(val,'cfg_dep')
    % Check dependency match
    sts2 = cellfun(@(cspec)match(item,cspec),{val.tgt_spec});
    if ~all(sts2)
        cfg_message('matlabbatch:checkval', ...
            '%s: Dependency does not match.', subsasgn_checkstr(item,subs));
    end
    val = val(sts2);
    sts = any(sts2);
else
    switch item.strtype
        case {'s'}
            if ~ischar(val)
                cfg_message('matlabbatch:checkval:strtype', ...
                        '%s: Item must be a string.', subsasgn_checkstr(item,substruct('.','val')));
                sts = false;
            else
                [sts, val] = numcheck(item,val);
                if sts && ~isempty(item.extras) && (ischar(item.extras) || iscellstr(item.extras))
                    pats = cellstr(item.extras);
                    mch = regexp(val, pats);
                    sts = any(~cellfun(@isempty, mch));
                    if ~sts
                        cfg_message('matlabbatch:checkval:strtype', ...
                            '%s: Item must match one of these patterns:\n%s', subsasgn_checkstr(item,substruct('.','val')), sprintf('%s\n', pats{:}));
                        sts = false;
                    end
                end
            end
        case {'s+'}
            cfg_message('matlabbatch:checkval:strtype', ...
                    '%s: FAILURE: Cant do s+ yet', subsasgn_checkstr(item,substruct('.','val')));
        case {'f'}
            % test whether val is a function handle or a name of an
            % existing function
            sts = subsasgn_check_funhandle(val);
            if ~sts
                cfg_message('matlabbatch:checkval:strtype', ...
                        '%s: Item must be a function handle or function name.', ...
                        subsasgn_checkstr(item,substruct('.','val')));
            end
        case {'n'}
            tol = 4*eps;
            sts = isempty(val) || (isnumeric(val) && all(val(~isnan(val(:))) >= 1) && ...
                                   all(abs(round(val(isfinite(val(:))))-val(isfinite(val(:)))) <= tol));
            if ~sts
                cfg_message('matlabbatch:checkval:strtype', ...
                        '%s: Item must be an array of natural numbers.', subsasgn_checkstr(item,substruct('.','val')));
                return;
            end
            [sts, val] = numcheck(item,val);
        case {'i'}
            tol = 4*eps;
            sts = isempty(val) || (isnumeric(val) && ...
                                   all(abs(round(val(isfinite(val(:))))-val(isfinite(val(:)))) <= tol));
            if ~sts
                cfg_message('matlabbatch:checkval:strtype', ...
                        '%s: Item must be an array of integers.', subsasgn_checkstr(item,substruct('.','val')));
                return;
            end
            [sts, val] = numcheck(item,val);
        case {'r'}
            sts = isempty(val) || (isnumeric(val) && all(isreal(val(:))));
            if ~sts
                cfg_message('matlabbatch:checkval:strtype', ...
                        '%s: Item must be an array of real numbers.', subsasgn_checkstr(item,substruct('.','val')));
                return;
            end
            [sts, val] = numcheck(item,val);
        case {'w'}
            tol = 4*eps;
            sts = isempty(val) || (isnumeric(val) && all(val(~isnan(val(:))) >= 0) && ...
                                   all(abs(round(val(isfinite(val(:))))-val(isfinite(val(:)))) <= tol));
            if ~sts
                cfg_message('matlabbatch:checkval:strtype', ...
                        '%s: Item must be an array of whole numbers.', subsasgn_checkstr(item,substruct('.','val')));
                return;
            end
            [sts, val] = numcheck(item,val);
        case {'e'}
            if ~isempty(item.extras) && subsasgn_check_funhandle(item.extras)
                [sts, val] = feval(item.extras, val, item.num);
            else
                [sts, val] = numcheck(item,val);
            end
        otherwise
            % only do size check for other strtypes
            [sts, val] = numcheck(item,val);
    end
end

function [sts, val] = numcheck(item,val)
% allow arbitrary size, if num field is empty
sts = true;
csz = size(val);
if ~isempty(item.num)
    if item.strtype == 's' && numel(item.num) == 2
        % interpret num field as [min max] # elements
        sts = item.num(1) <= numel(val) && numel(val) <= item.num(2);
        if ~sts
            cfg_message('matlabbatch:checkval:numcheck:mismatch', ...
                    '%s: Size mismatch (required [%s], present [%s]).', ...
                    subsasgn_checkstr(item,substruct('.','val')), num2str(item.num), num2str(csz));
        end
    else
        ind = item.num>0 & isfinite(item.num);
        if numel(csz) == 2
            % also try transpose for 2D arrays
            cszt = size(val');
        else
            cszt = csz;
        end
        if numel(item.num) ~= numel(csz)
            cfg_message('matlabbatch:checkval:numcheck:mismatch', ...
                    '%s: Dimension mismatch (required %d, present %d).', subsasgn_checkstr(item,substruct('.','val')), numel(item.num), numel(csz));
            sts = false;
            return;
        end
        if any(item.num(ind)-csz(ind))
            if any(item.num(ind)-cszt(ind))
                cfg_message('matlabbatch:checkval:numcheck:mismatch', ...
                        '%s: Size mismatch (required [%s], present [%s]).', ...
                        subsasgn_checkstr(item,substruct('.','val')), num2str(item.num), num2str(csz));
                sts = false;
                return
            else
                val = val';
                cfg_message('matlabbatch:checkval:numcheck:transposed', ...
                        '%s: Value transposed to match required size [%s].', ...
                        subsasgn_checkstr(item,substruct('.','val')), num2str(item.num));
            end
        end
    end
end