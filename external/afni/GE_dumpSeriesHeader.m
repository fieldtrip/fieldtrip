function outstr = GE_dumpSeriesHeader(series)
%
% outstr = GE_dumpSeriesHeader(series)
%
% Writes the series header to strimg the outstr
%
% Souheil Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

outstr = sprintf('se_suid: %s\n'        , char(series.se_suid));
outstr = strcat(outstr, sprintf('se_uniq: %d\n'        , series.se_uniq));
outstr = strcat(outstr, sprintf('se_diskid: %s\n'      , char(series.se_diskid)));
outstr = strcat(outstr, sprintf('se_exno: %d\n'        , series.se_exno));
outstr = strcat(outstr, sprintf('se_no: %d\n'          , series.se_no));
outstr = strcat(outstr, sprintf('se_datetime: %d\n'    , series.se_datetime));
outstr = strcat(outstr, sprintf('se_actual_dt: %d\n'   , series.se_actual_dt));
outstr = strcat(outstr, sprintf('se_desc: %s\n'        , cstring(series.se_desc)));
outstr = strcat(outstr, sprintf('pr_sysid: %s\n'       , cstring(series.pr_sysid)));
outstr = strcat(outstr, sprintf('pansysid: %s\n'       , cstring(series.pansysid)));
outstr = strcat(outstr, sprintf('se_typ: %d\n'         , series.se_typ));
outstr = strcat(outstr, sprintf('se_source: %d\n'      , series.se_source));
outstr = strcat(outstr, sprintf('se_plane: %d\n'       , series.se_plane));
outstr = strcat(outstr, sprintf('scan_type: %d\n'      , series.scan_type));
outstr = strcat(outstr, sprintf('position: %d\n'       , series.position));
outstr = strcat(outstr, sprintf('entry: %d\n'          , series.entry));
outstr = strcat(outstr, sprintf('anref: %s\n'          , char(series.anref)));
outstr = strcat(outstr, sprintf('lmhor: %f\n'          , series.lmhor));
outstr = strcat(outstr, sprintf('prtcl: %s\n'          , cstring(series.prtcl)));
outstr = strcat(outstr, sprintf('se_contrast: %d\n'    , series.se_contrast));
outstr = strcat(outstr, sprintf('start_ras: %s\n'      , char(series.start_ras)));
outstr = strcat(outstr, sprintf('start_loc: %f\n'      , series.start_loc));
outstr = strcat(outstr, sprintf('end_ras: %s\n'        , char(series.end_ras)));
outstr = strcat(outstr, sprintf('end_loc: %f\n'        , series.end_loc));
outstr = strcat(outstr, sprintf('se_pseq: %d\n'        , series.se_pseq));
outstr = strcat(outstr, sprintf('se_sortorder: %d\n'   , series.se_sortorder));
outstr = strcat(outstr, sprintf('se_lndmrkcnt: %d\n'   , series.se_lndmrkcnt));
outstr = strcat(outstr, sprintf('se_nacq: %d\n'        , series.se_nacq));
outstr = strcat(outstr, sprintf('xbasest: %d\n'        , series.xbasest));
outstr = strcat(outstr, sprintf('xbaseend: %d\n'       , series.xbaseend));
outstr = strcat(outstr, sprintf('xenhst: %d\n'         , series.xenhst));
outstr = strcat(outstr, sprintf('xenhend: %d\n'        , series.xenhend));
outstr = strcat(outstr, sprintf('se_lastmod: %d\n'     , series.se_lastmod));
outstr = strcat(outstr, sprintf('se_alloc_key: %s\n'   , cstring(series.se_alloc_key)));
outstr = strcat(outstr, sprintf('se_delta_cnt: %d\n'   , series.se_delta_cnt));
outstr = strcat(outstr, sprintf('se_verscre: %s\n'     , char(series.se_verscre)));
outstr = strcat(outstr, sprintf('se_verscur: %s\n'     , char(series.se_verscur)));
outstr = strcat(outstr, sprintf('se_pds_a: %f\n'       , series.se_pds_a));
outstr = strcat(outstr, sprintf('se_pds_c: %f\n'       , series.se_pds_c));
outstr = strcat(outstr, sprintf('se_pds_u: %f\n'       , series.se_pds_u));
outstr = strcat(outstr, sprintf('se_checksum: %d\n'    , series.se_checksum));
outstr = strcat(outstr, sprintf('se_complete: %d\n'    , series.se_complete));
outstr = strcat(outstr, sprintf('se_numarch: %d\n'     , series.se_numarch));
outstr = strcat(outstr, sprintf('se_imagect: %d\n'     , series.se_imagect));
outstr = strcat(outstr, sprintf('se_numimages: %d\n'   , series.se_numimages));
outstr = strcat(outstr, sprintf('se_images.length: %d\n'      , series.se_images.length));
outstr = strcat(outstr, sprintf('se_images.data: %d\n'      , series.se_images.data));
outstr = strcat(outstr, sprintf('se_numunimg: %d\n'    , series.se_numunimg));
outstr = strcat(outstr, sprintf('se_unimages.length: %d\n'    , series.se_unimages.length));
outstr = strcat(outstr, sprintf('se_unimages.data: %d\n'    , series.se_unimages.data));
outstr = strcat(outstr, sprintf('se_toarchcnt: %d\n'   , series.se_toarchcnt));
outstr = strcat(outstr, sprintf('se_toarchive.length: %d\n'   , series.se_toarchive.length));
outstr = strcat(outstr, sprintf('se_toarchive.data: %d\n'   , series.se_toarchive.data));
outstr = strcat(outstr, sprintf('echo1_alpha: %f\n'    , series.echo1_alpha));
outstr = strcat(outstr, sprintf('echo1_beta: %f\n'     , series.echo1_beta));
outstr = strcat(outstr, sprintf('echo1_window: %d\n'   , series.echo1_window));
outstr = strcat(outstr, sprintf('echo1_level: %d\n'    , series.echo1_level));
outstr = strcat(outstr, sprintf('echo2_alpha: %f\n'    , series.echo2_alpha));
outstr = strcat(outstr, sprintf('echo2_beta: %f\n'     , series.echo2_beta));
outstr = strcat(outstr, sprintf('echo2_window: %d\n'   , series.echo2_window));
outstr = strcat(outstr, sprintf('echo2_level: %d\n'    , series.echo2_level));
outstr = strcat(outstr, sprintf('echo3_alpha: %f\n'    , series.echo3_alpha));
outstr = strcat(outstr, sprintf('echo3_beta: %f\n'     , series.echo3_beta));
outstr = strcat(outstr, sprintf('echo3_window: %d\n'   , series.echo3_window));
outstr = strcat(outstr, sprintf('echo3_level: %d\n'    , series.echo3_level));
outstr = strcat(outstr, sprintf('echo4_alpha: %f\n'    , series.echo4_alpha));
outstr = strcat(outstr, sprintf('echo4_beta: %f\n'     , series.echo4_beta));
outstr = strcat(outstr, sprintf('echo4_window: %d\n'   , series.echo4_window));
outstr = strcat(outstr, sprintf('echo4_level: %d\n'    , series.echo4_level));
outstr = strcat(outstr, sprintf('echo5_alpha: %f\n'    , series.echo5_alpha));
outstr = strcat(outstr, sprintf('echo5_beta: %f\n'     , series.echo5_beta));
outstr = strcat(outstr, sprintf('echo5_window: %d\n'   , series.echo5_window));
outstr = strcat(outstr, sprintf('echo5_level: %d\n'    , series.echo5_level));
outstr = strcat(outstr, sprintf('echo6_alpha: %f\n'    , series.echo6_alpha));
outstr = strcat(outstr, sprintf('echo6_beta: %f\n'     , series.echo6_beta));
outstr = strcat(outstr, sprintf('echo6_window: %d\n'   , series.echo6_window));
outstr = strcat(outstr, sprintf('echo6_level: %d\n'    , series.echo6_level));
outstr = strcat(outstr, sprintf('echo7_alpha: %f\n'    , series.echo7_alpha));
outstr = strcat(outstr, sprintf('echo7_beta: %f\n'     , series.echo7_beta));
outstr = strcat(outstr, sprintf('echo7_window: %d\n'   , series.echo7_window));
outstr = strcat(outstr, sprintf('echo7_level: %d\n'    , series.echo7_level));
outstr = strcat(outstr, sprintf('echo8_alpha: %f\n'    , series.echo8_alpha));
outstr = strcat(outstr, sprintf('echo8_beta: %f\n'     , series.echo8_beta));
outstr = strcat(outstr, sprintf('echo8_window: %d\n'   , series.echo8_window));
outstr = strcat(outstr, sprintf('echo8_level: %d\n'    , series.echo8_level));
outstr = strcat(outstr, sprintf('series_uid: %s\n'     , sprintf('%s',series.series_uid)));
outstr = strcat(outstr, sprintf('landmark_uid: %s\n'   , sprintf('%s',series.landmark_uid)));
outstr = strcat(outstr, sprintf('equipmnt_uid: %s\n'   , sprintf('%s',series.equipmnt_uid)));
outstr = strcat(outstr, sprintf('se_padding: %s\n'     , cstring(series.se_padding)));

return
