function varargout = subsref(item, subs)

% function varargout = subsref(item, subs)
% subscript references we have to deal with are:
% one level
% item.(field)   - i.e. struct('type',{'.'} ,'subs',{field})
% item(idx)      - i.e. struct('type',{'()'},'subs',{idx})
% two levels
% item(idx).(field)
%
% to be dealt with elsewhere
% item.(field){fidx}
% three levels
% item(idx).(field){fidx}
% This function is identical for all classes derived from cfg_item, but it
% needs to be present in the class folder to access fields added by the
% derived class.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsref.m 4864 2012-08-27 13:57:31Z volkmar $

rev = '$Rev: 4864 $'; %#ok

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
switch subs(1).type,
    case {'.'},
        if numel(item) > 1
            cfg_message('matlabbatch:subsref:multiref', ...
                  'Field reference to multiple items not allowed for cfg_item classes.');
        end;
        switch subs(1).subs
            case local_mysubs_fields,
                val = item.(subs(1).subs);
            case par_fields,
                val = item.(par_class).(subs(1).subs);
            otherwise
                cfg_message('matlabbatch:subsref:unknownfield', ...
                      ['Reference to unknown field ''%s''.\nTo reference ' ...
                       'a field in the job structure, use a reference like ' ...
                       '''(x).%s'''], subs(1).subs, subs(1).subs);
        end;
        if numel(subs) > 1
            [varargout{1:nargout}] = cfg_callbuiltin('subsref',val,subs(2:end));
        else
            varargout = {val};
        end
%    case {'()','{}'},
%        val = subsref_job(item, subs, false);
    otherwise
        cfg_message('matlabbatch:subsref:unknowntype', ...
              'Unknown subsref type: ''%s''. This should not happen.', subs(1).type);
end
