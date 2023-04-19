function lst = fiff_rename_list(lst, ch_rename)

me = 'MNE:fiff_rename_list';

if nargin ~= 2
    error(me,'Incorrect number of arguments');
end
if isempty(lst)
    return;
end
if ~iscell(lst)
    error(me,'input must be a cell array:%s',lst);
end

for k = 1:length(lst)
    name = lst(k);
    if length(ch_rename)
        idx = find(strcmp(name, ch_rename{:, 1}));
        if length(idx)
            name = ch_rename{idx(1), 2};
        end
    end
    lst(k) = name;
end

return

end
