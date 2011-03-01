function output = read_micromed_trc(filename, begsample, endsample)

%--------------------------------------------------------------------------
% reads Micromed .TRC file into matlab, version Mariska
% input: filename
% output: datamatrix
%--------------------------------------------------------------------------

% ---------------- Opening File------------------
fid=fopen(filename,'rb');
if fid==-1
  error('Can''t open *.trc file')
end

%------------------reading patient & recording info----------
fseek(fid,64,-1);
header.surname=char(fread(fid,22,'char'))';
header.name=char(fread(fid,20,'char'))';

fseek(fid,128,-1);
day=fread(fid,1,'char');
if length(num2str(day))<2
  day=['0' num2str(day)];
else
  day=num2str(day);
end
month=fread(fid,1,'char');
switch month
  case 1
    month='JAN';
  case 2
    month='FEB';
  case 3
    month='MAR';
  case 4
    month='APR';
  case 5
    month='MAY';
  case 6
    month='JUN';
  case 7
    month='JUL';
  case 8
    month='AUG';
  case 9
    month='SEP';
  case 10
    month='OCT';
  case 11
    month='NOV';
  case 12
    month='DEC';
end
header.day=day;
header.month=month;
header.year=num2str(fread(fid,1,'char')+1900);

%------------------ Reading Header Info ---------
fseek(fid,175,-1);
header.Header_Type=fread(fid,1,'char');
if header.Header_Type ~= 4
  error('*.trc file is not Micromed System98 Header type 4')
end

fseek(fid,138,-1);
header.Data_Start_Offset=fread(fid,1,'uint32');
header.Num_Chan=fread(fid,1,'uint16');
header.Multiplexer=fread(fid,1,'uint16');
header.Rate_Min=fread(fid,1,'uint16');
header.Bytes=fread(fid,1,'uint16');

fseek(fid,176+8,-1);
header.Code_Area=fread(fid,1,'uint32');
header.Code_Area_Length=fread(fid,1,'uint32');

fseek(fid,192+8,-1);
header.Electrode_Area=fread(fid,1,'uint32');
header.Electrode_Area_Length=fread(fid,1,'uint32');

fseek(fid,400+8,-1);
header.Trigger_Area=fread(fid,1,'uint32');
header.Tigger_Area_Length=fread(fid,1,'uint32');

%----------------- Read Trace Data ----------

if nargin==1
  % determine the number of samples
  fseek(fid,header.Data_Start_Offset,-1);
  datbeg = ftell(fid);
  fseek(fid,0,1);
  datend = ftell(fid);
  header.Num_Samples = (datend-datbeg)/(header.Bytes*header.Num_Chan);
  if rem(header.Num_Samples, 1)~=0
    warning('rounding off the number of samples');
    header.Num_Samples = floor(header.Num_Samples);
  end
  % output the header
  output = header;
else
  % determine the selection of data to read
  if isempty(begsample)
    begsample = 1;
  end
  if isempty(endsample) || isinf(endsample)
    endsample = header.Num_Samples;
  end
  fseek(fid,header.Data_Start_Offset,-1);
  fseek(fid, header.Num_Chan*header.Bytes*(begsample-1), 0);
  switch header.Bytes
    case 1
      data = fread(fid, [header.Num_Chan endsample-begsample+1], 'uint8');
    case 2
      data = fread(fid, [header.Num_Chan endsample-begsample+1], 'uint16');
    case 4
      data = fread(fid, [header.Num_Chan endsample-begsample+1], 'uint32');
  end
  % output the data
  output = data-32768;
  % FIXME why is this value of -32768 subtracted?
  % FIXME some sort of calibration should be applied to get it into microvolt
end
fclose(fid);
