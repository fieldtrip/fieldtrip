function b = subsref(a,S)
% B = A.FIELDNAME; % SUBSREF - subscripted indexing
%
% See also GET

% For structure indexing 
if length(S)==1
    switch S.type
    case '.'
        % Ex. a.field
        b = get(a,S.subs);
        return;
    end;
elseif length(S)==2
     % For a.Channels(hw_ch_num) only
    if strcmpi(S(1).type,'.')&strcmpi(S(2).type,'()') 
        if strcmpi(S(1).subs,'Channels')   
            b = builtin('subsref',a,S);
            return;
        end;
    end;
end;
error('Incorrect parameter structure given to subsref');

        