function ch_rename = fiff_make_ch_rename(chs)

me='MNE:fiff_make_ch_rename';
if nargin ~= 1
    error(me,'Incorrect number of arguments');
end

ch_rename = {};
need_rename = false;
counts = {};
for k = 1:length(chs)
    if length(chs(k).ch_name) > 15
        need_rename = true;
    end
end
if ~need_rename
    return;
end

% count how many of each 15-char one we have
for k = 1:length(chs)
    short_name = chs(k).ch_name;
    short_name = short_name(1:min(length(short_name),15));
    if length(counts) == 0 || ~any(strcmp(counts(:, 1), short_name))
        counts = [counts; {short_name, 0, 0}];
    end
    idx = find(strcmp(counts(:, 1), short_name));
    counts{idx(1), 2} = counts{idx(1), 2} + 1;
end
% now do the assignments, taking into account duplicates and adding -1, -2, etc
for k = 1:length(chs)
    short_name = chs(k).ch_name;
    if length(short_name) > 15
        short_name = short_name(1:15);
        idx = find(strcmp(counts(:, 1), short_name));
        n = counts{idx(1), 2};
        if n > 1
            fw = ceil(log10(n));
            n = counts{idx(1), 3} + 1;
            counts{idx(1), 3} = n;
            short_name = sprintf(sprintf('%%s-%%0%dd', fw), short_name(1:15-fw-1), n);
        end
        ch_rename = [ch_rename; {chs(k).ch_name, short_name}];
    end
end

return;
end
