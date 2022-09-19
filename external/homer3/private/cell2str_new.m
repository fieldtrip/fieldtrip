function s = cell2str_new(c)
s = c;
if iscell(c)
    maxlen = 0;
    for ii = 1:length(c)
        if length(c{ii}) > maxlen
            maxlen = length(c{ii});
        end
    end
    s = char(zeros(length(c), maxlen));
    for ii = 1:length(c)
        s(ii,1:length(c{ii})) = c{ii};
    end
end

