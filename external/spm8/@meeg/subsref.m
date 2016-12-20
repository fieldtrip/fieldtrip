function varargout=subsref(this,subs)
% SUBSREF Subscripted reference
% An overloaded function...
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Stefan Kiebel
% $Id: subsref.m 2716 2009-02-09 17:14:15Z vladimir $

if isempty(subs)
    return;
end;

if this.Nsamples == 0
    error('Attempt to reference a field of an empty meeg object.');
end;

switch subs(1).type
    case '()'
        if numel(subs)~= 1, error('Expression too complicated');end;
        varargout = {double(subsref(this.data.y, subs))};
    case '{}'
    case '.'
        if ismethod(this, subs(1).subs)
            if numel(subs) == 1
                varargout = {feval(subs(1).subs, this)};
            elseif (numel(subs) == 2) && isequal(subs(2).type,  '()')
                varargout = {feval(subs(1).subs, this, subs(2).subs{:})};
            elseif (numel(subs)> 2) && isequal(subs(2).type,  '()')
                varargout{1} = builtin('subsref', ...
                    feval(subs(1).subs, this, subs(2).subs{:}),  subs(3:end));
            else
                varargout{1} = builtin('subsref', feval(subs(1).subs, this),  subs(2:end));
            end
        elseif isfield(this.other, subs(1).subs)
            field = getfield(this.other, subs(1).subs);
            if numel(subs)==1
                varargout = {field};
            else
                varargout{1} = builtin('subsref', field, subs(2:end));
            end
        else
            error('Reference to non-existent or private meeg method or field.');
        end
    otherwise
        error('Unfamiliar referencing type');
end


