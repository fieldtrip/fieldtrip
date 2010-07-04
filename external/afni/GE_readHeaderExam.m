function ex_hdr = GE_readHeaderExam(fid, byte_align)
%
% ex_hdr = GE_readHeaderExam(fid, byte_align)
%
% Loads the exam header from the file with filed id fid
% and returns it as a structure. 
% if byte_align = 1 then 32-bit alignment (SGI, LX2 format)
% if byte_align = 0 then 16-bit alignment (Sun, 5.X format)
%
%
% Souheil J. Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%


% define the structure and read ifseek(fid,1,0);  % byte alignmentn the data
% to overcome the byte alignment problems
% break up the assignment into pieces using the setfield function
ex_hdr = struct('ex_suid', fread(fid,4,'uchar')); %Suite ID for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_uniq', fread(fid,1,'int16'));        %The Make-Unique Flag%
ex_hdr = setfield(ex_hdr, 'ex_diskid', fread(fid,1,'uchar'));      %Disk ID for this Exam%
fseek(fid,1,0); % 16-bit alignment
ex_hdr = setfield(ex_hdr, 'ex_no', fread(fid,1,'uint16'));         %Exam Number%
ex_hdr = setfield(ex_hdr, 'hospname', fread(fid,33,'uchar'));      %Hospital Name%
fseek(fid,1,0); % 16-bit alignment
ex_hdr = setfield(ex_hdr, 'detect', fread(fid,1,'int16'));         %Detector Type%
if byte_align; fseek(fid,2,0); end % 32-bit alignment
ex_hdr = setfield(ex_hdr, 'numcells', fread(fid,1,'int32'));       %Number of cells in det%
ex_hdr = setfield(ex_hdr, 'zerocell', fread(fid,1,'float32'));     %Cell number at theta%
ex_hdr = setfield(ex_hdr, 'cellspace', fread(fid,1,'float32'));    %Cell spacing%
ex_hdr = setfield(ex_hdr, 'srctodet', fread(fid,1,'float32'));     %Distance from source to detector%
ex_hdr = setfield(ex_hdr, 'srctoiso', fread(fid,1,'float32'));     %Distance from source to iso%
ex_hdr = setfield(ex_hdr, 'tubetyp', fread(fid,1,'int16'));        %Tube type%
ex_hdr = setfield(ex_hdr, 'dastyp', fread(fid,1,'int16'));         %DAS type%
ex_hdr = setfield(ex_hdr, 'num_dcnk', fread(fid,1,'int16'));       %Number of Decon Kernals%
ex_hdr = setfield(ex_hdr, 'dcn_len', fread(fid,1,'int16'));        %Number of elements in a Decon Kernal%
ex_hdr = setfield(ex_hdr, 'dcn_density', fread(fid,1,'int16'));    %Decon Kernal density%
ex_hdr = setfield(ex_hdr, 'dcn_stepsize', fread(fid,1,'int16'));   %Decon Kernal stepsize%
ex_hdr = setfield(ex_hdr, 'dcn_shiftcnt', fread(fid,1,'int16'));   %Decon Kernal Shift Count%
if byte_align; fseek(fid,2,0); end % 32-bit alignment
ex_hdr = setfield(ex_hdr, 'magstrength', fread(fid,1,'int32'));    %Magnet strength (in gauss)%
ex_hdr = setfield(ex_hdr, 'patid', fread(fid,13,'uchar'));         %Patient ID for this Exam%
ex_hdr = setfield(ex_hdr, 'patname', fread(fid,25,'uchar'));       %Patientsda Name%
ex_hdr = setfield(ex_hdr, 'patage', fread(fid,1,'int16'));         %Patient Age (years, months or days)%
ex_hdr = setfield(ex_hdr, 'patian', fread(fid,1,'int16'));         %Patient Age Notation%
ex_hdr = setfield(ex_hdr, 'patsex', fread(fid,1,'int16'));         %Patient Sex%
ex_hdr = setfield(ex_hdr, 'patweight', fread(fid,1,'int32'));      %Patient Weight%
ex_hdr = setfield(ex_hdr, 'trauma', fread(fid,1,'int16'));         %Trauma Flag%
ex_hdr = setfield(ex_hdr, 'hist', fread(fid,61,'uchar'));          %Patient History%
ex_hdr = setfield(ex_hdr, 'reqnum', fread(fid,13,'uchar'));        %Requisition Number%
ex_hdr = setfield(ex_hdr, 'ex_datetime', fread(fid,1,'int32'));    %Exam date/time stamp%
ex_hdr = setfield(ex_hdr, 'refphy', fread(fid,33,'uchar'));        %Referring Physician%
ex_hdr = setfield(ex_hdr, 'diagrad', fread(fid,33,'uchar'));       %Diagnostician/Radiologist%
ex_hdr = setfield(ex_hdr, 'op', fread(fid,4,'uchar'));             %Operator%
ex_hdr = setfield(ex_hdr, 'ex_desc', fread(fid,23,'uchar'));       %Exam Description%
ex_hdr = setfield(ex_hdr, 'ex_typ', fread(fid,3,'uchar'));         %Exam Type%
ex_hdr = setfield(ex_hdr, 'ex_format', fread(fid,1,'int16'));      %Exam Format%
if byte_align; fseek(fid,6,0); end % 32-bit alignment
ex_hdr = setfield(ex_hdr, 'firstaxtime', fread(fid,1,'float64'));  %Start time(secs) of first axial in exam%
ex_hdr = setfield(ex_hdr, 'ex_sysid', fread(fid,9,'uchar'));       %Creator Suite and Host%
fseek(fid,1,0); % 16-bit alignment
if byte_align; fseek(fid,2,0); end % 32-bit alignment
ex_hdr = setfield(ex_hdr, 'ex_lastmod', fread(fid,1,'int32'));     %Date/Time of Last Change%
ex_hdr = setfield(ex_hdr, 'protocolflag', fread(fid,1,'int16'));   %Non-Zero indicates Protocol Exam%
ex_hdr = setfield(ex_hdr, 'ex_alloc_key', fread(fid,13,'uchar'));  %Process that allocated this record%
fseek(fid,1,0); % 16-bit alignment
ex_hdr = setfield(ex_hdr, 'ex_delta_cnt', fread(fid,1,'int32'));   %Indicates number of updates to header%
ex_hdr = setfield(ex_hdr, 'ex_verscre', fread(fid,2,'uchar'));     %Genesis Version - Created%
ex_hdr = setfield(ex_hdr, 'ex_verscur', fread(fid,2,'uchar'));     %Genesis Version - Now%
ex_hdr = setfield(ex_hdr, 'ex_checksum', fread(fid,1,'uint32'));   %Exam Record Checksum%
ex_hdr = setfield(ex_hdr, 'ex_complete', fread(fid,1,'int32'));    %Exam Complete Flag%
ex_hdr = setfield(ex_hdr, 'ex_seriesct', fread(fid,1,'int32'));    %Last Series Number Used%
ex_hdr = setfield(ex_hdr, 'ex_numarch', fread(fid,1,'int32'));     %Number of Series Archived%
ex_hdr = setfield(ex_hdr, 'ex_numseries', fread(fid,1,'int32'));   %Number of Series Existing%
ex_hdr = setfield(ex_hdr, 'ex_series', struct('length', fread(fid,1,'uint32'), ...
                                              'data', fread(fid,1,'uint32')));  %Series Keys for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_numunser', fread(fid,1,'int32'));    %Number of Unstored Series%
ex_hdr = setfield(ex_hdr, 'ex_unseries', struct('length', fread(fid,1,'uint32'), ...
                                                'data', fread(fid,1,'uint32')));    %Unstored Series Keys for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_toarchcnt', fread(fid,1,'int32'));   %Number of Unarchived Series%
ex_hdr = setfield(ex_hdr, 'ex_toarchive', struct('length', fread(fid,1,'uint32'), ...
                                                 'data', fread(fid,1,'uint32')));     %Unarchived Series Keys for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_prospcnt', fread(fid,1,'int32'));    %Number of Prospective/Scout Series%
ex_hdr = setfield(ex_hdr, 'ex_prosp', struct('length', fread(fid,1,'uint32'), ...
                                             'data', fread(fid,1,'uint32'))); %Prospective/Scout Series Keys for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_modelnum', fread(fid,1,'int32'));    %Last Model Number used%
ex_hdr = setfield(ex_hdr, 'ex_modelcnt', fread(fid,1,'int32'));    %Number of ThreeD Models%
ex_hdr = setfield(ex_hdr, 'ex_models', struct('length', fread(fid,1,'uint32'), ...
                                              'data', fread(fid,1,'uint32')));  %ThreeD Model Keys for Exam%
ex_hdr = setfield(ex_hdr, 'ex_stat', fread(fid,1,'int16'));        %Patient Status%
ex_hdr = setfield(ex_hdr, 'uniq_sys_id', fread(fid,16,'uchar'));   %Unique System ID%
ex_hdr = setfield(ex_hdr, 'service_id', fread(fid,16,'uchar'));    %Unique Service ID%
ex_hdr = setfield(ex_hdr, 'mobile_loc', fread(fid,4,'uchar'));     %Mobile Location Number%
ex_hdr = setfield(ex_hdr, 'study_uid', fread(fid,32,'uchar'));     %Study Entity Unique ID%
ex_hdr = setfield(ex_hdr, 'study_status', fread(fid,1,'int16'));   %indicates if study has complete info(DICOM/genesis)%
ex_hdr = setfield(ex_hdr, 'ex_padding', fread(fid,516,'uchar'));   %Spare Space%
if byte_align; fseek(fid,4,0); end % byte alignment

return
