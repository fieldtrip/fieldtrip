function json=moxunit_fieldtrip_util_struct2json(s)
% convert struct to json string representation
%
% json=moxunit_fieldtrip_util_struct2json(s)
%
% Inputs:
%   s                   struct for which the values must all be strings,
%                       or scalar logical or double arrays
%
% Output:
%   json                JSON string representation of the input
%

    if ~isstruct(s)
        error('input must be a struct');
    end

    keys=fieldnames(s);
    n=numel(keys);
    args=cell(2,n);
    for k=1:n
        key=keys{k};
        raw_value=s.(key);

        if ischar(raw_value)
            value=sprintf('"%s"',escape(raw_value));

        elseif isnumeric(raw_value) && numel(raw_value)==1
            value=sprintf('%d',raw_value);

        elseif islogical(raw_value) && numel(raw_value)==1
            if raw_value
                value='true';
            else
                value='false';
            end

        else
            error('Unsupported input for key %s: <class %s>',...
                            key, class(raw_value));
        end

        args{1,k}=sprintf('"%s":%s',key,value);

        if k<n
            sep=',';
        else
            sep='';
        end
        args{2,k}=sep;
    end

    json=sprintf('{%s}',sprintf('%s',args{:}));

function s=escape(s)
    s=regexprep(s,'"','\\"');
