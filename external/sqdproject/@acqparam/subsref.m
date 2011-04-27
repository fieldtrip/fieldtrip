function b = subsref(a,S)
% B = A.FIELDNAME; % SUBSREF - subscripted indexing
%
% SEE ALSO, GET

% For structure indexing 
switch S.type
% Ex. a.field
case '.'
    b = get(a,S.subs);
end;

        