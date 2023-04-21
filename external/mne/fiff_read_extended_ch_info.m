function [chs, ch_rename] = fiff_read_extended_ch_info(chs,meas_info,fid)

me='MNE:fiff_read_extended_ch_info';
if nargin ~= 3
    error(me,'Incorrect number of arguments');
end

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

ch_rename = {};

ch_infos = fiff_dir_tree_find(meas_info,FIFF.FIFFB_CH_INFO);

if length(ch_infos) == 0
    return;
elseif length(ch_infos) ~= length(chs)
    error('extended channel info does not match channels')
end
for k=1:length(chs)
    new = ch_infos(k);
    ch = chs(k);
    for p = 1:new.nent
        kind = new.dir(p).kind;
        if kind == FIFF.FIFF_CH_DACQ_NAME
            tag = fiff_read_tag(fid, new.dir(p).pos);
            data = tag.data;
            ch_rename = [ch_rename; {ch.ch_name, data}];
            ch.ch_name = data;
            break
        end
    end
    chs(k) = ch;
end

return;

end
