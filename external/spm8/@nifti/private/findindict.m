function entry = findindict(c,dcode)
% Look up an entry in the dictionary
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: findindict.m 1143 2008-02-07 19:33:33Z spm $


entry = [];
d = getdict;
d = d.(dcode);
if ischar(c)
    for i=1:length(d),
        if strcmpi(d(i).label,c),
            entry = d(i);
            break;
        end;
    end;
elseif isnumeric(c) && numel(c)==1
    for i=1:length(d),
        if d(i).code==c,
            entry = d(i);
            break;
        end;
    end;
else
    error(['Inappropriate code for ' dcode '.']);
end;
if isempty(entry)
    fprintf('\nThis is not an option.  Try one of these:\n');
    for i=1:length(d)
        fprintf('%5d) %s\n', d(i).code, d(i).label);
    end;
    %fprintf('\nNO CHANGES MADE\n');
end;

