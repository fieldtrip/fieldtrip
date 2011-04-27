function a = subsasgn(a,S,val)
% A.FIELD = val; %Subscripted assignment
%
% SEE ALSO, SET

% For structure array
switch S.type
case '.'
    if nargin == 3 % Ex. a.field = val;
        a = set(a,S.subs,val);
    end;
end;
        