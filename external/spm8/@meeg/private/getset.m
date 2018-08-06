function res = getset(this, parent, fieldname, ind, values)
% Generic method for getting and setting multiple fields of meeg struct
% FORMAT res = getset(this, parent, fieldname, ind, values)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: getset.m 4430 2011-08-12 18:47:17Z vladimir $

this = struct(this);

if nargin == 3 || isempty(ind) || (ischar(ind) && isequal(ind, ':'))
    try
        ind = 1:numel(getfield(this, parent));
    catch
        res = [];
        return;
    end
end

% Get
if nargin <= 4
    res = cell(1, length(ind));
    for i = 1:length(ind)
        res{i} = getfield(this, parent, {ind(i)}, fieldname);
    end

    if all(cellfun('isclass', res, 'double') & ~cellfun('isempty', res))
        res = [res{:}];
    end

    if iscell(res) && (numel(res) == 1) && (numel(getfield(this, parent)) == 1) &&...
            strcmp(parent, 'trials') && strcmp(this.type, 'continuous')
        res = res{1};
    end
    
    return
end

% Set
if nargin == 5          
    % This might fail in some pathological cases, but not in what it's
    % supposed to be used for.
    if (isnumeric(values) || islogical(values)) && (length(values) == length(ind))
        values = num2cell(values);
    end
    
    if iscell(values) && ~(numel(values)==1 || numel(values) == length(ind))
        error('Illegal assignment: cannot match values and indices.');
    end        
        
    for i = 1:length(ind)
        if iscell(values)
            this = setfield(this, parent, {ind(i)}, fieldname, values{i});
        else
            this = setfield(this, parent, {ind(i)}, fieldname, values);
        end
    end
    res = meeg(this);
end