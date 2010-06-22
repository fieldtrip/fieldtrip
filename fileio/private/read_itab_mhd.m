function mhd = read_itab_mhd(filename)

fid = fopen(filename, 'rb');

%  Name of structure
mhd.stname = fread(fid, [1 10], 'uint8=>char');           %  Header identifier  (VP_BIOMAG)
mhd.stver = fread(fid, [1 8], 'uint8=>char');             %  Header version
mhd.stendian = fread(fid, [1 4], 'uint8=>char');          %  LE (little endian) or BE (big endian) format

%  Subject's INFOs
mhd.first_name = fread(fid, [1 32], 'uint8=>char');       %  Subject's first name
mhd.last_name = fread(fid, [1 32], 'uint8=>char');        %  Subject's family name
mhd.id = fread(fid, [1 32], 'uint8=>char');               %  Subject's id
mhd.notes = fread(fid, [1 256], 'uint8=>char');           %  Notes on measurement

%  Other subj infos
mhd.subj_info.sex         = fread(fid, [1 1], 'uint8=>char');  %  Sex (M or F)
pad = fread(fid, [1 5], 'uint8');
mhd.subj_info.notes       = fread(fid, [1 256], 'uint8=>char'); %  Notes on subject
mhd.subj_info.height      = fread(fid, [1 1], 'float');  %  Height in cm
mhd.subj_info.weight      = fread(fid, [1 1], 'float');  %  Weight in kg
mhd.subj_info.birthday    = fread(fid, [1 1], 'int32'); %  Birtday (1-31)
mhd.subj_info.birthmonth  = fread(fid, [1 1], 'int32'); %  Birthmonth (1-12)
mhd.subj_info.birthyear   = fread(fid, [1 1], 'int32');  %  Birthyear (1900-2002)

%  Data acquisition INFOs
mhd.time = fread(fid, [1 12], 'uint8=>char');             %  time (ascii)
mhd.date = fread(fid, [1 16], 'uint8=>char');             %  date (ascii)

mhd.nchan     = fread(fid, [1 1], 'int32');          %  total number of channels
mhd.nelech    = fread(fid, [1 1], 'int32');          %  number of electric channels
mhd.nelerefch = fread(fid, [1 1], 'int32');          %  number of electric reference channels
mhd.nmagch    = fread(fid, [1 1], 'int32');          %  number of magnetic channels
mhd.nmagrefch = fread(fid, [1 1], 'int32');          %  number of magnetic reference channels
mhd.nauxch    = fread(fid, [1 1], 'int32');          %  number of auxiliary channels
mhd.nparamch  = fread(fid, [1 1], 'int32');          %  number of parameter channels
mhd.ndigitch  = fread(fid, [1 1], 'int32');          %  number of digit channels
mhd.nflagch   = fread(fid, [1 1], 'int32');          %  number of flag channels

mhd.data_type   = fread(fid, [1 1], 'int32');        %  0 - BE_SHORT      (HP-PA, big endian)
%  1 - BE_LONG       (HP-PA, big endian)
%  2 - BE_FLOAT      (HP-PA, big endian)
%  3 - LE_SHORT      (Intel, little endian)
%  4 - LE_LONG       (Intel, little endian)
%  5 - LE_FLOAT      (Intel, little endian)
%  6 - RTE_A_SHORT   (HP-A900, big endian)
%  7 - RTE_A_FLOAT   (HP-A900, big endian)
%  8 - ASCII

mhd.smpfq       = fread(fid, [1 1], 'single');        %  sampling frequency in Hz
mhd.hw_low_fr   = fread(fid, [1 1], 'single');        %  hw data acquisition low pass filter
mhd.hw_hig_fr   = fread(fid, [1 1], 'single');        %  hw data acquisition high pass filter
mhd.hw_comb     = fread(fid, [1 1], 'int32');        %  hw data acquisition 50 Hz filter  (1-TRUE)
mhd.sw_hig_tc   = fread(fid, [1 1], 'single');        %  sw data acquisition high pass time constant
mhd.compensation= fread(fid, [1 1], 'int32');        %  0 - no compensation    1 - compensation
mhd.ntpdata     = fread(fid, [1 1], 'int32');        %  total number of time points in data

mhd.no_segments    = fread(fid, [1 1], 'int32');     %  Number of segments described in the segment structure

%  INFOs on different segments
for i=1:5
  mhd.sgmt(i).start = fread(fid, [1 1], 'int32');       %  Starting time point from beginning of data
  mhd.sgmt(i).ntptot = fread(fid, [1 1], 'int32');      %  Total number of time points
  mhd.sgmt(i).type = fread(fid, [1 1], 'int32');
  mhd.sgmt(i).no_samples = fread(fid, [1 1], 'int32');
  mhd.sgmt(i).st_sample = fread(fid, [1 1], 'int32');
end

mhd.nsmpl          = fread(fid, [1 1], 'int32');     %  Overall number of samples

%  INFOs on different samples
for i=1:4096
  mhd.smpl(i).start = fread(fid, [1 1], 'int32');       %  Starting time point from beginning of data
  mhd.smpl(i).ntptot = fread(fid, [1 1], 'int32');      %  Total number of time points
  mhd.smpl(i).ntppre = fread(fid, [1 1], 'int32');      %  Number of points in pretrigger
  mhd.smpl(i).type = fread(fid, [1 1], 'int32');
  mhd.smpl(i).quality = fread(fid, [1 1], 'int32');
end


mhd.nrefchan        = fread(fid, [1 1], 'int32');    %  number of reference channels
mhd.ref_ch = fread(fid, [1 640], 'int32');          %  reference channel list
mhd.ntpref          = fread(fid, [1 1], 'int32');    %  total number of time points in reference

%  Header INFOs
mhd.raw_header_type   = fread(fid, [1 1], 'int32');  %  0  - Unknown header
%  2  - rawfile (A900)
%  3  - GE runfile header
%  31 - ATB runfile header version 1.0
%  41 - IFN runfile header version 1.0
%  51 - BMDSys runfile header version 1.0

mhd.header_type     = fread(fid, [1 1], 'int32');    %  0  - Unknown header
%  2  - rawfile (A900)
%  3  - GE runfile header
%  4  - old header
%  10 - 256ch normal header
%  11 - 256ch master header
%  20 - 640ch normal header
%  21 - 640ch master header
%  31 - ATB runfile header version 1.0
%  41 - IFN runfile header version 1.0
%  51 - BMDSys runfile header version 1.0

mhd.conf_file = fread(fid, [1 64], 'uint8=>char');        %  Filename used for data acquisition configuration

mhd.header_size      = fread(fid, [1 1], 'int32');   %  sizeof(header) at the time of file creation
mhd.start_reference  = fread(fid, [1 1], 'int32');   %  start reference
mhd.start_data       = fread(fid, [1 1], 'int32');   %  start data

mhd.rawfile          = fread(fid, [1 1], 'int32');   %  0 - not a rawfile      1 - rawfile
mhd.multiplexed_data = fread(fid, [1 1], 'int32');   %  0 - FALSE              1 - TRUE

mhd.isns        = fread(fid, [1 1], 'int32');     %  sensor code   1 - Single channel
%               28 - Original Rome 28 ch.
%               29 -  ..............
%               .. -  ..............
%               45 - Updated Rome 28 ch. (spring 2009)
%               .. -  ..............
%               .. -  ..............
%               55 - Original Chieti 55 ch. flat
%              153 - Original Chieti 153 ch. helmet
%              154 - Chieti 153 ch. helmet from Jan 2002

%  Channel's INFOs
for i=1:640
  mhd.ch(i).type = fread(fid, [1 1], 'uint8');
  pad = fread(fid, [1 3], 'uint8');
  %  type      0 - unknown
  %            1 - ele
  %            2 - mag
  %            4 - ele ref
  %            8 - mag ref
  %           16 - aux
  %           32 - param
  %           64 - digit
  %          128 - flag
  mhd.ch(i).number = fread(fid, [1 1], 'int32');         %  number
  mhd.ch(i).label = fixstr(fread(fid, [1 16], 'uint8=>char'));      %  label
  mhd.ch(i).flag = fread(fid, [1 1], 'uint8');
  pad = fread(fid, [1 3], 'uint8');
  %  on/off flag 0 - working channel
  %              1 - noisy channel
  %              2 - very noisy channel
  %              3 - broken channel
  mhd.ch(i).amvbit = fread(fid, [1 1], 'float');         %  calibration from LSB to mV
  mhd.ch(i).calib = fread(fid, [1 1], 'float');          %  calibration from mV to unit
  mhd.ch(i).unit = fread(fid, [1 6], 'uint8=>char');  %  unit label  (fT, uV, ...)
  pad = fread(fid, [1 2], 'uint8');
  mhd.ch(i).ncoils = fread(fid, [1 1], 'int32');         %  number of coils building up one channel
  mhd.ch(i).wgt = fread(fid, [1 10], 'float');       %  weight of coils
  %  position and orientation of coils
  for j=1:10
    mhd.ch(i).position(j).r_s = fread(fid, [1 3], 'float');
    mhd.ch(i).position(j).u_s = fread(fid, [1 3], 'float');
  end
end


%  Sensor position INFOs

mhd.r_center = fread(fid, [1 3], 'float');          %  sensor position in convenient format
mhd.u_center = fread(fid, [1 3], 'float');
mhd.th_center      = fread(fid, [1 1], 'float');  %  sensor orientation as from markers fit
mhd.fi_center= fread(fid, [1 1], 'float');
mhd.rotation_angle= fread(fid, [1 1], 'float');
mhd.cosdir = fread(fid, [3 3], 'float');         %  for compatibility only
mhd.irefsys      = fread(fid, [1 1], 'int32');    %  reference system   0 - sensor reference system
%                     1 - Polhemus
%                     2 - head3
%                     3 - MEG

%  Marker positions for MRI integration
mhd.num_markers  = fread(fid, [1 1], 'int32');  %  Total number of markers
mhd.i_coil = fread(fid, [1 64], 'int32');      %  Markers to be used to find sensor position
mhd.marker = fread(fid, [3 64], 'int32');      %  Position of all the markers
mhd.best_chi     = fread(fid, [1 1], 'float');  %  Best chi_square value obtained in finding sensor position

mhd.cup_vertex_center  = fread(fid, [1 1], 'float');    %  dist anc sensor  cente vertex  (as entered
%  from keyboard)
mhd.cup_fi_center  = fread(fid, [1 1], 'float');        %  fi angle of sensor center
mhd.cup_rotation_angle  = fread(fid, [1 1], 'float');   %  rotation angle of sensor center axis

mhd.dist_a1_a2  = fread(fid, [1 1], 'float');           %  head informations
%  (used to find subject's head dimensions)
mhd.dist_inion_nasion = fread(fid, [1 1], 'float');
mhd.max_circ = fread(fid, [1 1], 'float');
mhd.nasion_vertex_inion = fread(fid, [1 1], 'float');

%  Data analysis INFOs
mhd.security  = fread(fid, [1 1], 'int32');          %  security flag
mhd.ave_alignement  = fread(fid, [1 1], 'int32');    %  average data alignement 0 - FALSE
%                          1 - TRUE

mhd.itri  = fread(fid, [1 1], 'int32');              %  trigger channel number
mhd.ntpch  = fread(fid, [1 1], 'int32');             %  no. of time points per channel
mhd.ntppre  = fread(fid, [1 1], 'int32');            %  no. of time points of pretrigger
mhd.navrg  = fread(fid, [1 1], 'int32');             %  no. of averages
mhd.nover  = fread(fid, [1 1], 'int32');             %  no. of discarded averages
mhd.nave_filt  = fread(fid, [1 1], 'int32');         %  no. of applied filters

%  Filters used before average
for i=1:15
  mhd.ave_filt(i).type = fread(fid, [1 1], 'int32');                 % type   0 - no filter
  %        1 - bandpass  - param[0]: highpass freq
  %                      - param[1]: lowpass freq
  %        2 - notch     - param[0]: notch freq 1
  %                        param[1]: notch freq 2
  %                        param[2]: notch freq 3
  %                        param[3]: notch freq 4
  %                        param[4]: span
  %        3 - artifact  - param[0]: True/False
  %        4 - adaptive  - param[0]: True/False
  %        5 - rectifier - param[0]: True/False
  %        6 - heart     - param[0]: True/False
  %        7 - evoked    - param[0]: True/False
  %        8 - derivate  - param[0]: True/False
  %        9 - polarity  - param[0]: True/False
  mhd.ave_filt(i).param = fread(fid, [1 5], 'float');             % up to 5 filter parameters
end

mhd.stdev  = fread(fid, [1 1], 'int32');             %  0 - not present
mhd.bas_start  = fread(fid, [1 1], 'int32');         %  starting data points for baseline
mhd.bas_end  = fread(fid, [1 1], 'int32');           %  ending data points for baseline

mhd.source_files  = fread(fid, [1 32], 'int32');  %  Progressive number of files (if more than one)

%  template INFOs
mhd.ichtpl = fread(fid, [1 1], 'int32');
mhd.ntptpl = fread(fid, [1 1], 'int32');
mhd.ifitpl = fread(fid, [1 1], 'int32');
mhd.corlim = fread(fid, [1 1], 'float');

%  Filters used before template
for i=1:15
  mhd.tpl_filt(i).type = fread(fid, [1 1], 'int32');                 % type   0 - no filter
  %        1 - bandpass  - param[0]: highpass freq
  %                      - param[1]: lowpass freq
  %        2 - notch     - param[0]: notch freq 1
  %                        param[1]: notch freq 2
  %                        param[2]: notch freq 3
  %                        param[3]: notch freq 4
  %                        param[4]: span
  %        3 - artifact  - param[0]: True/False
  %        4 - adaptive  - param[0]: True/False
  %        5 - rectifier - param[0]: True/False
  %        6 - heart     - param[0]: True/False
  %        7 - evoked    - param[0]: True/False
  %        8 - derivate  - param[0]: True/False
  %        9 - polarity  - param[0]: True/False
  mhd.tpl_filt(i).param = fread(fid, [1 5], 'float');             % up to 5 filter parameters
end

%  Just in case info
mhd.dummy = fread(fid, [1 64], 'int32');
% there seems to be more dummy data at the end...

fclose(fid);

function str = fixstr(str)
sel = find(str==0, 1, 'first');
if ~isempty(sel)
  str = str(1:sel-1);
end

