function a = subsasgn(a,S,val)
% A.FIELD = val; %Subscripted assignment
%
% SEE ALSO, SET

if isempty(a)
    % Ex. a(1) = chanhandle(filename,channel_num);
    if isequal(S.type, '()')&isa(val,'chanhandle')&isequal(S.subs{:},1)
        a = val;
        return;
    else
        newline= sprintf('\n');
        str = ['Could not initialize chanhandle object with given parameters',newline,...
        'Please note that:',newline,...
        '  1. Right-hand side should be a valid chanhandle object',newline,...
        '  2. Only the first channel should be set',newline,...
        '  3. Property of chanhandle cannot be set until the object itself is initialzed'];
        error(str);        
    end;
end;
if nargin == 3
    switch length(S)
    case 1
        switch S.type
        case '()'
            if isa(val,'chanhandle')
                if isequal(S.subs{1}, 1:length(val)) % a(i) = chanhandle_obj;
                    a.handle(S.subs{1}) = getfield(struct(val),'handle');
                elseif isequal(length(S.subs{:}),length(getfield(struct(val),'handle')))
                    % Ex. a(i1,i2,...) = [chanhandle_obj1,chanhandle_obj2,...];
                    for i = 1:length(S.subs)
                        a.handle(S.subs{i}) = getfield(struct(val(i)),'handle');
                    end;
                else
                    error('Input chanhandle structure and Channel number do not match');
                end;
            elseif isstruct(val)  %Ex. a(i) = chanhandle_obj.handle;
                if isequal(length(S.subs{:}),length(val))
                    for i = 1:length(S.subs)
                        a.handle(S.subs{i}) = val(i);
                    end;
                else
                    error('Input chanhandle structure and Channel number do not match');
                end;
            else
                error('RHS must be either a valid structure or a chanhandle object');
            end;                    
        case '.'
            a = set(a,S.subs,val);      % a.field = val;
        end;
    case 2
        if strcmp(S(1).type,'()')&strcmpi(S(2).type,'.')    % a(i).field = val;
            a = set(a,S(1).subs{:},S(2).subs,val);
        elseif strcmp(S(1).type,'.')&strcmpi(S(2).type,'()')% a.field(i) = val;
            a = set(a,S(1).subs,S(2).subs{:},val);
        elseif strcmp(S(1).type,'.')&strcmpi(S(2).type,'.') % a.field1.field2 = val;
            a = set(a,[S(1).subs,'.',S(2).subs],val);
        else
            error('Sorry - output structure not supported');
        end;
    otherwise
        error('Sorry - output structure not supported');
    end;
else
    error('incorrect number of arguments');
end;