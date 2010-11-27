function se_hdr = GE_readHeaderSeries(fid, byte_align)
%
% se_hdr = read_series_header(fid, byte_align)
%
% Loads the series header from a file with filed id fid
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

% define the structure and read in the data
% to overcome the byte alignment problems
% break up the assignment into pieces using the setfield function
se_hdr = struct('se_suid', fread(fid,4,'uchar'));                       %Suite ID for this Series%
se_hdr = setfield(se_hdr, 'se_uniq', fread(fid,1,'int16'));            %The Make-Unique Flag%
se_hdr = setfield(se_hdr, 'se_diskid', fread(fid,1,'uchar'));          %Disk ID for this Series%
fseek(fid,1,0); % 16-bit alignment
se_hdr = setfield(se_hdr, 'se_exno', fread(fid,1,'uint16'));            %Exam Number%
se_hdr = setfield(se_hdr, 'se_no ', fread(fid,1,'int16'));              %Series Number%
se_hdr = setfield(se_hdr, 'se_datetime', fread(fid,1,'int32'));        %Allocation Series Data/Time stamp%
se_hdr = setfield(se_hdr, 'se_actual_dt', fread(fid,1,'int32'));       %Actual Series Data/Time stamp%
se_hdr = setfield(se_hdr, 'se_desc', fread(fid,30,'uchar'));       %Series Description%
se_hdr = setfield(se_hdr, 'pr_sysid', fread(fid,9,'uchar'));       %Primary Receiver Suite and Host%
se_hdr = setfield(se_hdr, 'pansysid', fread(fid,9,'uchar'));       %Archiver Suite and Host%
se_hdr = setfield(se_hdr, 'se_typ ', fread(fid,1,'int16'));             %Series Type%
se_hdr = setfield(se_hdr, 'se_source ', fread(fid,1,'int16'));          %Series from which prescribed%
se_hdr = setfield(se_hdr, 'se_plane ', fread(fid,1,'int16'));           %Most-like Plane (for L/S)%
se_hdr = setfield(se_hdr, 'scan_type ', fread(fid,1,'int16'));          %Scout or Axial (for CT)%
se_hdr = setfield(se_hdr, 'position', fread(fid,1,'int32'));           %Patient Position%
se_hdr = setfield(se_hdr, 'entry', fread(fid,1,'int32'));              %Patient Entry%
se_hdr = setfield(se_hdr, 'anref', fread(fid,3,'uchar'));          %Anatomical reference%
fseek(fid,1,0); % 16-bit alignment
se_hdr = setfield(se_hdr, 'lmhor', fread(fid,1,'float32'));              %Horizontal Landmark%
se_hdr = setfield(se_hdr, 'prtcl', fread(fid,25,'uchar'));         %Scan Protocol Name%
fseek(fid,1,0); % 16-bit alignment
se_hdr = setfield(se_hdr, 'se_contrast ', fread(fid,1,'int16'));        %Non-zero if > 0 image used contrast(L/S)%
se_hdr = setfield(se_hdr, 'start_ras', fread(fid,1,'uchar'));          %RAS letter for first scan location (L/S)%
fseek(fid,1,0); % 16-bit alignment
if byte_align; fseek(fid,2,0); end % 32-bit alignment
se_hdr = setfield(se_hdr, 'start_loc', fread(fid,1,'float32'));          %First scan location (L/S)%
se_hdr = setfield(se_hdr, 'end_ras', fread(fid,1,'uchar'));            %RAS letter for last scan location (L/S)%
fseek(fid,1,0); % 16-bit alignment
if byte_align; fseek(fid,2,0); end % 32-bit alignment
se_hdr = setfield(se_hdr, 'end_loc', fread(fid,1,'float32'));            %Last scan location (L/S)%
se_hdr = setfield(se_hdr, 'se_pseq ', fread(fid,1,'int16'));            %Last Pulse Sequence Used (L/S)%
se_hdr = setfield(se_hdr, 'se_sortorder ', fread(fid,1,'int16'));       %Image Sort Order (L/S)%
se_hdr = setfield(se_hdr, 'se_lndmrkcnt', fread(fid,1,'int32'));       %Landmark Counter%
se_hdr = setfield(se_hdr, 'se_nacq ', fread(fid,1,'int16'));            %Number of Acquisitions%
se_hdr = setfield(se_hdr, 'xbasest ', fread(fid,1,'int16'));            %Starting number for baselines%
se_hdr = setfield(se_hdr, 'xbaseend', fread(fid,1,'int16'));            %Ending number for baselines%
se_hdr = setfield(se_hdr, 'xenhst', fread(fid,1,'int16'));             %Starting number for enhanced scans%
se_hdr = setfield(se_hdr, 'xenhend', fread(fid,1,'int16'));            %Ending number for enhanced scans%
if byte_align; fseek(fid,2,0); end % 32-bit alignment
se_hdr = setfield(se_hdr, 'se_lastmod', fread(fid,1,'int32'));         %Date/Time of Last Change%
se_hdr = setfield(se_hdr, 'se_alloc_key', fread(fid,13,'uchar'));  %Process that allocated this record%
fseek(fid,1,0); % 16-bit alignment
if byte_align; fseek(fid,2,0); end % 32-bit alignment
se_hdr = setfield(se_hdr, 'se_delta_cnt', fread(fid,1,'int32'));       %Indicates number of updates to header%
se_hdr = setfield(se_hdr, 'se_verscre', fread(fid,2,'uchar'));     %Genesis Version - Created%
se_hdr = setfield(se_hdr, 'se_verscur', fread(fid,2,'uchar'));     %Genesis Version - Now%
se_hdr = setfield(se_hdr, 'se_pds_a', fread(fid,1,'float32'));           %PixelData size - as stored%
se_hdr = setfield(se_hdr, 'se_pds_c', fread(fid,1,'float32'));           %PixelData size - Compressed%
se_hdr = setfield(se_hdr, 'se_pds_u', fread(fid,1,'float32'));           %PixelData size - UnCompressed%
se_hdr = setfield(se_hdr, 'se_checksum', fread(fid,1,'uint32'));        %Series Record checksum%
se_hdr = setfield(se_hdr, 'se_complete', fread(fid,1,'int32'));        %Series Complete Flag%
se_hdr = setfield(se_hdr, 'se_numarch', fread(fid,1,'int32'));         %Number of Images Archived%
se_hdr = setfield(se_hdr, 'se_imagect', fread(fid,1,'int32'));         %Last Image Number Used%
se_hdr = setfield(se_hdr, 'se_numimages', fread(fid,1,'int32'));       %Number of Images Existing%
se_hdr = setfield(se_hdr, 'se_images', struct('length', fread(fid,1,'uint32'), ...
                                               'data', fread(fid,1,'uint32'))); %Image Keys for this Series%
se_hdr = setfield(se_hdr, 'se_numunimg', fread(fid,1,'int32'));        %Number of Unstored Images%
se_hdr = setfield(se_hdr, 'se_unimages', struct('length', fread(fid,1,'uint32'), ...
                                               'data', fread(fid,1,'uint32'))); %Unstored Image Keys for this Series%
se_hdr = setfield(se_hdr, 'se_toarchcnt', fread(fid,1,'int32'));       %Number of Unarchived Images%
se_hdr = setfield(se_hdr, 'se_toarchive', struct('length', fread(fid,1,'uint32'), ...
                                               'data', fread(fid,1,'uint32'))); %Unarchived Image Keys for this Series%
se_hdr = setfield(se_hdr, 'echo1_alpha', fread(fid,1,'float32'));        %Echo 1 Alpha Value%
se_hdr = setfield(se_hdr, 'echo1_beta', fread(fid,1,'float32'));         %Echo 1 Beta Value%
se_hdr = setfield(se_hdr, 'echo1_window', fread(fid,1,'uint16'));       %Echo 1 Window Value%
se_hdr = setfield(se_hdr, 'echo1_level', fread(fid,1,'int16'));        %Echo 1 Level Value%
se_hdr = setfield(se_hdr, 'echo2_alpha', fread(fid,1,'float32'));        %Echo 2 Alpha Value%
se_hdr = setfield(se_hdr, 'echo2_beta', fread(fid,1,'float32'));         %Echo 2 Beta Value%
se_hdr = setfield(se_hdr, 'echo2_window', fread(fid,1,'uint16'));       %Echo 2 Window Value%
se_hdr = setfield(se_hdr, 'echo2_level', fread(fid,1,'int16'));        %Echo 2 Level Value%
se_hdr = setfield(se_hdr, 'echo3_alpha', fread(fid,1,'float32'));        %Echo 3 Alpha Value%
se_hdr = setfield(se_hdr, 'echo3_beta', fread(fid,1,'float32'));         %Echo 3 Beta Value%
se_hdr = setfield(se_hdr, 'echo3_window', fread(fid,1,'uint16'));       %Echo 3 Window Value%
se_hdr = setfield(se_hdr, 'echo3_level', fread(fid,1,'int16'));        %Echo 3 Level Value%
se_hdr = setfield(se_hdr, 'echo4_alpha', fread(fid,1,'float32'));        %Echo 4 Alpha Value%
se_hdr = setfield(se_hdr, 'echo4_beta', fread(fid,1,'float32'));         %Echo 4 Beta Value%
se_hdr = setfield(se_hdr, 'echo4_window', fread(fid,1,'uint16'));       %Echo 4 Window Value%
se_hdr = setfield(se_hdr, 'echo4_level', fread(fid,1,'int16'));        %Echo 4 Level Value%
se_hdr = setfield(se_hdr, 'echo5_alpha', fread(fid,1,'float32'));        %Echo 5 Alpha Value%
se_hdr = setfield(se_hdr, 'echo5_beta', fread(fid,1,'float32'));         %Echo 5 Beta Value%
se_hdr = setfield(se_hdr, 'echo5_window', fread(fid,1,'uint16'));       %Echo 5 Window Value%
se_hdr = setfield(se_hdr, 'echo5_level', fread(fid,1,'int16'));        %Echo 5 Level Value%
se_hdr = setfield(se_hdr, 'echo6_alpha', fread(fid,1,'float32'));        %Echo 6 Alpha Value%
se_hdr = setfield(se_hdr, 'echo6_beta', fread(fid,1,'float32'));         %Echo 6 Beta Value%
se_hdr = setfield(se_hdr, 'echo6_window', fread(fid,1,'uint16'));       %Echo 6 Window Value%
se_hdr = setfield(se_hdr, 'echo6_level', fread(fid,1,'int16'));        %Echo 6 Level Value%
se_hdr = setfield(se_hdr, 'echo7_alpha', fread(fid,1,'float32'));        %Echo 7 Alpha Value%
se_hdr = setfield(se_hdr, 'echo7_beta', fread(fid,1,'float32'));         %Echo 7 Beta Value%
se_hdr = setfield(se_hdr, 'echo7_window', fread(fid,1,'uint16'));       %Echo 7 Window Value%
se_hdr = setfield(se_hdr, 'echo7_level', fread(fid,1,'int16'));        %Echo 7 Level Value%
se_hdr = setfield(se_hdr, 'echo8_alpha', fread(fid,1,'float32'));        %Echo 8 Alpha Value%
se_hdr = setfield(se_hdr, 'echo8_beta', fread(fid,1,'float32'));         %Echo 8 Beta Value%
se_hdr = setfield(se_hdr, 'echo8_window', fread(fid,1,'uint16'));       %Echo 8 Window Value%
se_hdr = setfield(se_hdr, 'echo8_level', fread(fid,1,'int16'));        %Echo 8 Level Value%
se_hdr = setfield(se_hdr, 'series_uid', fread(fid,32,'uchar'));    %Series Entity Unique ID%
se_hdr = setfield(se_hdr, 'landmark_uid', fread(fid,32,'uchar'));  %Landmark Unique ID%
se_hdr = setfield(se_hdr, 'equipmnt_uid', fread(fid,32,'uchar'));  %Equipment Unique ID%
se_hdr = setfield(se_hdr, 'se_padding', fread(fid,588,'uchar'));   %Spare Space%

return
