function cals=fiff_write_ch_infos(fid,chs,reset_range,ch_rename)

me='MNE:fiff_write_ch_infos';

if nargin == 3
    ch_rename = {};
elseif nargin ~= 4
    error(me,'Incorrect number of arguments');
end

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

write_rename = false;
for k = 1:length(chs)
    %
    %   Scan numbers may have been messed up
    %
    ch = chs(k);
    ch.scanno = k;
    if reset_range
        ch.range = 1.0;
    end
    cals(k) = ch.cal;

    if length(ch_rename)
        idx = find(strcmp(ch_rename{:, 1}, ch.ch_name));
        if length(idx)
            write_rename = true;
            ch.ch_name = ch_rename{idx(1), 2};
        end
    end
    fiff_write_ch_info(fid,ch);
end

% add extra struct
if write_rename
    for k=1:length(chs)
        ch = chs(k);
        fiff_start_block(fid,FIFF.FIFFB_CH_INFO);
        fiff_write_string(fid,FIFF.FIFF_CH_DACQ_NAME,ch.ch_name);
        fiff_write_int(fid, FIFF.FIFF_CH_SCAN_NO, ch.scanno);
        fiff_write_int(fid, FIFF.FIFF_CH_LOGICAL_NO, ch.logno);
        fiff_write_int(fid,FIFF.FIFF_CH_KIND,ch.kind);
        fiff_write_float(fid,FIFF.FIFF_CH_RANGE,ch.range);
        fiff_write_float(fid,FIFF.FIFF_CH_CAL,ch.cal);
        fiff_write_int(fid,FIFF.FIFF_CH_COIL_TYPE,ch.coil_type);
        fiff_write_float(fid,FIFF.FIFF_CH_LOC,ch.loc);
        fiff_write_int(fid,FIFF.FIFF_CH_UNIT,ch.unit);
        fiff_write_int(fid,FIFF.FIFF_CH_UNIT_MUL,ch.unit_mul);
        fiff_write_string(fid,FIFF.FIFF_CH_DACQ_NAME, ch.ch_name);
        fiff_write_int(fid,FIFF.FIFF_CH_COORD_FRAME,ch.coord_frame);
        fiff_end_block(fid, FIFF.FIFFB_CH_INFO)
    end
end

return;

end
