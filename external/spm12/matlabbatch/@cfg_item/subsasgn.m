function item = subsasgn(item, subs, varargin)

% function item = subsasgn(item, subs, varargin)
% This function implements subsasgn for all classes derived from cfg_item.
% It relies on the capability of each class constructor to re-classify a
% struct object after a new value has been assigned to its underlying
% struct (This capability has to be implemented in the derived class).
% The structure of a configuration tree does not permit any arrays of
% cfg_item objects. Therefore, the only subscript reference and
% assignment within an cfg_item is a dot assignment to fields of this
% cfg_item. 
% Subscript references we have to deal with are:
% one level
% item.(field)   - i.e. struct('type',{'.'} ,'subs',{field})
%
% to be dealt with elsewhere
% item.(field){fidx}
% 
% In a future version, '()' and '{}' subscripts may be supported to
% access val fields of a cfg_item tree as if they were part of a
% harvested job. For cfg_branch objects (where dot assignments are used
% for val fields in their job tree) it is mandatory to index the job as a
% struct array to access harvested fields.
% This function is identical for all classes derived from cfg_item. A
% copy of it must be present in each derived class to be able to access
% derived fields.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

persistent local_mysubs_fields;
persistent par_class;
persistent par_fields;
if ~iscell(local_mysubs_fields)
    local_mysubs_fields = mysubs_fields;
    citem = class(item);
    switch citem
        case 'cfg_exbranch',
            par_class = 'cfg_branch';
            pf1 = subs_fields(item.cfg_branch);
            pf2 = subs_fields(cfg_item);
            par_fields = [pf1(:); pf2(:)]';
        case 'cfg_item',
            par_class = '';
            par_fields = {};
        otherwise
            par_class = 'cfg_item';
            par_fields = subs_fields(item.cfg_item);
    end;
end
if numel(item) ~= 1
    cfg_message('matlabbatch:subsasgn', ...
          'Arrays of cfg_item objects not supported.');
end;

switch subs(1).type,
    case {'.'},
        if numel(item) == 1 && nargin == 3
            if numel(subs) > 1
                % Only part of field value(s) assigned. Get old field value(s),
                % assign values to it.
                if any(strcmp(subs(1).subs,par_fields))
                    val0 = item.(par_class).(subs(1).subs);
                else
                    val0 = item.(subs(1).subs);
                end
                val1 = cfg_callbuiltin('subsasgn',val0,subs(2:end),varargin{1});
            else
                % Field values to assign
                val1 = varargin{1};
            end
            [ok, val] = subsasgn_check(item,subs,val1);
            if ok,
                switch subs(1).subs
                    case par_fields,
                        item.(par_class).(subs(1).subs) = val;
                    case local_mysubs_fields,
                        item.(subs(1).subs) = val;
                    otherwise
                        cfg_message('matlabbatch:subsasgn', ...
                                    ['Reference to unknown field ''%s''.\n' ...
                                     'To assign to a field in the job structure, use a reference like ' ...
                                     '''(x).%s''.'], subs(1).subs, subs(1).subs);
                end;
            end;
        else
            cfg_message('matlabbatch:subsasgn', ...
                        ['Array assignments not supported for fields of cfg_item objects.\n' ...
                         'To assign to a field in the job structure, use a reference like ' ...
                         '''(x).%s''.'], subs(1).subs);

        end;
    case {'()','{}'},
        cfg_message('matlabbatch:subsasgn', ...
                    'Subscript type ''%s'' reserved for future use.', subs(1).type);
    otherwise
        cfg_message('matlabbatch:subsasgn', ...
                    'Unknown subsref type: ''%s''. This should not happen.', subs(1).type);
end;
return;