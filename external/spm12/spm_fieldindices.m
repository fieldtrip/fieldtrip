function [ix] = spm_fieldindices(X,varargin)
% Returns the indices of fields in a structure (and vice versa)
% FORMAT [i]     = spm_fieldindices(X,field1,field2,...)
% FORMAT [field] = spm_fieldindices(X,i)
%
% X         - structure
% field1,.. - fields
%
% i         - vector of indices or fieldname{s}
%
% Note: Fields are returned in column order of X, regardless of the order
% of fields specified in the input.
%
%__________________________________________________________________________
% Copyright (C) 2010-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fieldindices.m 6809 2016-06-15 12:51:55Z peter $


% if varargin is a vector simply return fieldnames
%--------------------------------------------------------------------------
if isnumeric(varargin{1})
    if numel(varargin{1}) > 1
        for j = 1:length(varargin{1})
            ix{j} = spm_fieldindices(X,varargin{1}(j));
        end
        return
    end
end

% create structure of zeros
%--------------------------------------------------------------------------
X0    = spm_zeros(X);
ix    = spm_vec(X0);

% and add one to specified fields
%--------------------------------------------------------------------------
for i = 1:length(varargin)
    
    if ischar(varargin{i})
        
        field = varargin{i};
        x  = X0;
        try
            % field is a field name
            %--------------------------------------------------------------
            f  = x.(field);
            f  = spm_unvec(spm_vec(f) + 1,f);
            x.(field) = f;
            ix = ix + spm_vec(x);
            
        catch
            try
                % field is an expression
                %----------------------------------------------------------
                eval(['x.' field ' = x.' field ' + 1;']);
                ix = ix + spm_vec(x);
            end
        end
        
    % or return the name of the field
    %----------------------------------------------------------------------
    elseif isnumeric(varargin{i})
        
        ind    = varargin{i};
        x      = ix;
        x(ind) = 1;
        x      = spm_unvec(x,X);
        name   = fieldnames(x);
        
        for j = 1:length(name)
            if any(spm_vec(x.(name{j})))
                if iscell(x.(name{j}))
                    for k = 1:numel(x.(name{j}))
                        if any(spm_vec(x.(name{j}){k}))
                            [p,q] = find(x.(name{j}){k});
                            ix = sprintf('%s{%i}(%i,%i)',name{j},k,p,q);
                            return
                        end
                    end
                elseif isnumeric(x.(name{j}))
                    s   = size(x.(name{j}));
                    if numel(s) < 3
                        if min(s) == 1
                            p = find(x.(name{j}));
                            ix = sprintf('%s(%i)',name{j},p);
                            return
                        else
                            [p,q] = find(x.(name{j}));
                            ix = sprintf('%s(%i,%i)',name{j},p,q);
                            return
                        end
                    else
                        [p,q,r] = ind2sub(s,find(x.(name{j})));
                        ix = sprintf('%s(%i,%i,%i)',name{j},p,q,r);
                        return
                    end
                end
            end
        end
    end
end

% find indices
%--------------------------------------------------------------------------
ix = find(ix);
