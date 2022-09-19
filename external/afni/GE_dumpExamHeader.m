function outstr = GE_dumpExamHeader(exam)
%
% outstr = GE_dumpExamHeader(exam)
%
% Writes the exam header to the string outstr
%
% Souheil J. Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

outstr = sprintf('ex_suid: %s\n'     , char(exam.ex_suid));
outstr = strcat(outstr, sprintf('ex_uniq: %d\n'     , exam.ex_uniq));
outstr = strcat(outstr, sprintf('ex_diskid: %s\n'   , char(exam.ex_diskid)));
outstr = strcat(outstr, sprintf('ex_no: %d\n'       , exam.ex_no));
outstr = strcat(outstr, sprintf('hospname: %s\n'    , cstring(exam.hospname)));
outstr = strcat(outstr, sprintf('detect: %d\n'      , exam.detect));
outstr = strcat(outstr, sprintf('numcells: %d\n'    , exam.numcells));
outstr = strcat(outstr, sprintf('zerocell: %f\n'    , exam.zerocell));
outstr = strcat(outstr, sprintf('cellspace: %f\n'   , exam.cellspace));
outstr = strcat(outstr, sprintf('srctodet: %f\n'    , exam.srctodet));
outstr = strcat(outstr, sprintf('srctoiso: %f\n'    , exam.srctoiso));
outstr = strcat(outstr, sprintf('tubetyp: %d\n'     , exam.tubetyp));
outstr = strcat(outstr, sprintf('dastyp: %d\n'      , exam.dastyp));
outstr = strcat(outstr, sprintf('num_dcnk: %d\n'    , exam.num_dcnk));
outstr = strcat(outstr, sprintf('dcn_len: %d\n'     , exam.dcn_len));
outstr = strcat(outstr, sprintf('dcn_density: %d\n' , exam.dcn_density));
outstr = strcat(outstr, sprintf('dcn_stepsize: %d\n', exam.dcn_stepsize));
outstr = strcat(outstr, sprintf('dcn_shiftcnt: %d\n', exam.dcn_shiftcnt));
outstr = strcat(outstr, sprintf('magstrength: %d\n' , exam.magstrength));
outstr = strcat(outstr, sprintf('patid: %s\n'       , cstring(exam.patid)));
outstr = strcat(outstr, sprintf('patname: %s\n'     , cstring(exam.patname)));
outstr = strcat(outstr, sprintf('patage: %d\n'      , exam.patage));
outstr = strcat(outstr, sprintf('patian: %d\n'      , exam.patian));
outstr = strcat(outstr, sprintf('patsex: %d\n'      , exam.patsex));
outstr = strcat(outstr, sprintf('patweight: %d\n'   , exam.patweight));
outstr = strcat(outstr, sprintf('trauma: %d\n'      , exam.trauma));
outstr = strcat(outstr, sprintf('hist: %s\n'        , cstring(exam.hist)));
outstr = strcat(outstr, sprintf('reqnum: %s\n'      , cstring(exam.reqnum)));
outstr = strcat(outstr, sprintf('ex_datetime: %d\n' , exam.ex_datetime));
outstr = strcat(outstr, sprintf('refphy: %s\n'      , cstring(exam.refphy)));
outstr = strcat(outstr, sprintf('diagrad: %s\n'     , cstring(exam.diagrad)));
outstr = strcat(outstr, sprintf('op: %s\n'          , cstring(exam.op)));
outstr = strcat(outstr, sprintf('ex_desc: %s\n'     , cstring(exam.ex_desc)));
outstr = strcat(outstr, sprintf('ex_typ: %s\n'      , cstring(exam.ex_typ)));
outstr = strcat(outstr, sprintf('ex_format: %d\n'   , exam.ex_format));
outstr = strcat(outstr, sprintf('firstaxtime: %f\n' , exam.firstaxtime));
outstr = strcat(outstr, sprintf('ex_sysid: %s\n'    , sprintf('%s',exam.ex_sysid)));
outstr = strcat(outstr, sprintf('ex_lastmod: %d\n'  , exam.ex_lastmod));
outstr = strcat(outstr, sprintf('protocolflag: %d\n', exam.protocolflag));
outstr = strcat(outstr, sprintf('ex_alloc_key: %s\n', sprintf('%s',exam.ex_alloc_key)));
outstr = strcat(outstr, sprintf('ex_delta_cnt: %d\n', exam.ex_delta_cnt));
outstr = strcat(outstr, sprintf('ex_verscre: %s\n'  , char(exam.ex_verscre)));
outstr = strcat(outstr, sprintf('ex_verscur: %s\n'  , char(exam.ex_verscur)));
outstr = strcat(outstr, sprintf('ex_checksum: %d\n' , exam.ex_checksum));
outstr = strcat(outstr, sprintf('ex_complete: %d\n' , exam.ex_complete));
outstr = strcat(outstr, sprintf('ex_seriesct: %d\n' , exam.ex_seriesct));
outstr = strcat(outstr, sprintf('ex_numarch: %d\n'  , exam.ex_numarch));
outstr = strcat(outstr, sprintf('ex_numseries: %d\n', exam.ex_numseries));
outstr = strcat(outstr, sprintf('ex_series.length: %d\n'   , exam.ex_series.length));
outstr = strcat(outstr, sprintf('ex_series.data: %d\n'   , exam.ex_series.data));
outstr = strcat(outstr, sprintf('ex_numunser: %d\n' , exam.ex_numunser));
outstr = strcat(outstr, sprintf('ex_unseries.length: %d\n' , exam.ex_unseries.length));
outstr = strcat(outstr, sprintf('ex_unseries.data: %d\n' , exam.ex_unseries.data));
outstr = strcat(outstr, sprintf('ex_toarchcnt: %d\n', exam.ex_toarchcnt));
outstr = strcat(outstr, sprintf('ex_toarchive.length: %d\n', exam.ex_toarchive.length));
outstr = strcat(outstr, sprintf('ex_toarchive.data: %d\n', exam.ex_toarchive.data));
outstr = strcat(outstr, sprintf('ex_prospcnt: %d\n' , exam.ex_prospcnt));
outstr = strcat(outstr, sprintf('ex_prosp.length: %d\n'    , exam.ex_prosp.length));
outstr = strcat(outstr, sprintf('ex_prosp.data: %d\n'    , exam.ex_prosp.data));
outstr = strcat(outstr, sprintf('ex_modelnum: %d\n' , exam.ex_modelnum));
outstr = strcat(outstr, sprintf('ex_modelcnt: %d\n' , exam.ex_modelcnt));
outstr = strcat(outstr, sprintf('ex_models.length: %d\n'   , exam.ex_models.length));
outstr = strcat(outstr, sprintf('ex_models.data: %d\n'   , exam.ex_models.data));
outstr = strcat(outstr, sprintf('ex_stat: %d\n'     , exam.ex_stat));
outstr = strcat(outstr, sprintf('uniq_sys_id: %s\n' , sprintf('%s',exam.uniq_sys_id)));
outstr = strcat(outstr, sprintf('service_id: %s\n'  , sprintf('%s',exam.service_id)));
outstr = strcat(outstr, sprintf('mobile_loc: %s\n'  , sprintf('%s',exam.mobile_loc)));
outstr = strcat(outstr, sprintf('study_uid: %s\n'   , sprintf('%s',exam.study_uid)));
outstr = strcat(outstr, sprintf('study_status: %d\n', exam.study_status));
outstr = strcat(outstr, sprintf('ex_padding: %s\n' , cstring(exam.ex_padding)));

return
