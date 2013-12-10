function varargout = subsref(this,subs)
% Subscript referencing for GIfTI objects
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id$

if length(this) > 1 && ~strcmp(subs(1).type,'()')
    warning('Not implemented.');
    for i=1:numel(this)
        varargout{i} = subsref(this(i),subs);
    end
    return;
end

switch subs(1).type
    case '.'
        [i,j] = isintent(this,subs(1).subs);
        if isempty(i)
            if strcmp(subs(1).subs,'private')
                varargout{1} = builtin('struct',this);
            else
                error('Reference to non-existent field ''%s''.',subs(1).subs);
            end
        else
            if strcmp(subs(1).subs,'mat')
                varargout{1} = this.data{j}.space.MatrixData;
            elseif strcmp(subs(1).subs,'labels')
                varargout{1} = this.label;
            else
                if length(j) == 1
                    varargout{1} = this.data{j}.data;
                else
                    a = [this.data{j}];
                    try
                        varargout{1} = [a.data];
                    catch
                        error('Data arrays are of different sizes.');
                    end
                end
            end
        end
        if strcmp(subs(1).subs,'faces')
            varargout{1} = varargout{1} + 1; % indices start at 1
        end
        if length(subs) > 1
            varargout{1} = subsref(varargout{1},subs(2:end));
        end
    case '{}'
        error('Cell contents reference from a non-cell array object.');
    case '()'
        if length(subs) == 1
            varargout{1} = builtin('subsref',this,subs(1));
        else
            varargout{1} = subsref(builtin('subsref',this,subs(1)),subs(2:end));
        end
    otherwise
        error('This should not happen.');
end