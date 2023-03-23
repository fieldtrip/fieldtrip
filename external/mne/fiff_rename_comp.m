function comp = fiff_rename_comp(comp, ch_rename)

me='MNE:fiff_rename_comp';

if nargin ~= 2
    error(me,'Incorrect number of arguments');
end

comp.data.row_names = fiff_rename_list(comp.data.row_names, ch_rename);
comp.data.col_names = fiff_rename_list(comp.data.col_names, ch_rename);

return

end
