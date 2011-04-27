function a = subsasgn(a,S,val)
% A.FIELD = val; %Subscripted assignment
%
% SEE ALSO, SET

% For structure array
switch S.type
case '.'
    % Ex. a.field = val
    if nargin = 3
        a = set(a,S.subs,val);
    else
        error('not enough arguments given to subsasgn');
    end;
end;
