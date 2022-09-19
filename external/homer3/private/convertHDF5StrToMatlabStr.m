function snew = convertHDF5StrToMatlabStr(s, type)
snew = '';

if ~exist('type', 'var') || isempty(type)
    type = 'char';
end

% Convert muti-row char array to cell string array
if size(s,1)==1 && strcmp(type,'char')
    snew = strtrim_improve(s);
else
    snew = cell(size(s,1),1);
    for ii = 1:size(s,1)
        snew{ii} = strtrim_improve(s(ii,:));
    end
end
