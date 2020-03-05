function sts = match(item, spec)

% sts = match(item, spec)
% Spec must be a cell array of struct arrays with one or more fields. Each
% struct must contain two fields - 'name' and 'value'.
% An item matches, if it has a field with the specified field name and the
% contents of this field equals the contents of spec.value. If the field
% name is 'class', an item matches, if its class name is equal to
% spec.value.
% Matches within each struct array are OR-concatenated, while matches
% between struct arrays are AND-concatenated.
% An empty spec always matches.
% Special matching rules for cfg_entries apply to the .strtype field.
% An item.strtype
% 'e' - matches any strtype
% 'n' - matches strtype 'n'
% 'w' - matches strtype 'n', 'w'
% 'i' - matches strtype 'n', 'w', 'i'
% 'r' - matches strtype 'n', 'w', 'i', 'r'
% Any other strtype matches only on equality.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: match.m 2675 2009-01-30 14:58:07Z volkmar $

rev = '$Rev: 2675 $'; %#ok

% match an empty spec
sts = true;

for k = 1:numel(spec)
    % Assume no match
    sts = false;
    for l = 1:numel(spec{k})
        switch spec{k}(l).name,
            case 'strtype',
                switch item.strtype
                    case 'e'
                        % always match strtype 'e'
                        sts = true;
                    case 'w'
                        sts = any(strcmp(spec{k}(l).value,{'n','w'}));
                    case 'i'
                        sts = any(strcmp(spec{k}(l).value,{'n','w','i'}));
                    case 'r'
                        sts = any(strcmp(spec{k}(l).value,{'n','w','i','r'}));
                    otherwise
                        % exact match
                        sts = strcmp(item.strtype, spec{k}(l).value);
                end
            case 'num',
                sts = true;
            case 'class'
                sts = strcmpi(spec{k}(l).value,class(item));
            otherwise
                spec1{1}(1) = spec{k}(l);
                sts = match(item.cfg_item, spec1);
        end;
        if sts
            % OR: success on first match
            break;
        end;
    end;
    if ~sts
        % AND: fail on first non-match
        break;
    end;
end;
