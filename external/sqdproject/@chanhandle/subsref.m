function b = subsref(a,S)
% B = A.FIELDNAME; % SUBSREF - subscripted indexing
%
% See also GET

switch length(S)
case 1
    switch S.type
    case '()'
        b = get(a,S.subs{:});
    case '.'
        b = get(a,S.subs);
    end;
case 2
    if strcmp(S(1).type,'()')&strcmpi(S(2).type,'.')    % a(i).field
        b = get(a,S(1).subs{:},S(2).subs);
    elseif strcmp(S(1).type,'.')&strcmpi(S(2).type,'()')% a.field(i)
        b = get(a,S(1).subs,S(2).subs{:});
    elseif strcmp(S(1).type,'.')&strcmpi(S(2).type,'.') % a.field1.field2
        b = get(a,[S(1).subs,'.',S(2).subs]);
    else
        error('Sorry - output structure not supported');
    end;
otherwise
    error('Sorry - output structure not supported');
end;