function output = read_micromed_trc(filename, begsample, endsample)

%--------------------------------------------------------------------------
% reads Micromed .TRC file into matlab, version Mariska, edited by Romain
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


% ============== Retrieving electrode info  ===============
% Order
fseek(fid,184,'bof');
OrderOff = fread(fid,1,'ulong');
fseek(fid,OrderOff,'bof');
vOrder = zeros(header.Num_Chan,1);
for iChan = 1 : header.Num_Chan, vOrder(iChan) = fread(fid,1,'ushort'); end
fseek(fid,200,'bof');
ElecOff = fread(fid,1,'ulong');
for iChan = 1 : header.Num_Chan
    fseek(fid,ElecOff+128*vOrder(iChan),'bof');
    if ~fread(fid,1,'uchar'), continue; end
    header.elec(iChan).bip = fread(fid,1,'uchar');
    header.elec(iChan).Name = deblank(char(fread(fid,6,'char'))');
    header.elec(iChan).Name(isspace(header.elec(iChan).Name)) = []; % remove spaces
    header.elec(iChan).Ref = deblank(char(fread(fid,6,'char'))');
    header.elec(iChan).LogicMin = fread(fid,1,'long');
    header.elec(iChan).LogicMax = fread(fid,1,'long');
    header.elec(iChan).LogicGnd = fread(fid,1,'long');
    header.elec(iChan).PhysMin = fread(fid,1,'long');
    header.elec(iChan).PhysMax = fread(fid,1,'long');
    Unit = fread(fid,1,'ushort');
    switch Unit
        case -1
            header.elec(iChan).Unit = 'nV';
        case 0
            header.elec(iChan).Unit = 'uV';
        case 1
            header.elec(iChan).Unit = 'mV';
        case 2
            header.elec(iChan).Unit = 'V';
        case 100
            header.elec(iChan).Unit = '%';
        case 101
            header.elec(iChan).Unit = 'bpm';
        case 102
            header.elec(iChan).Unit = 'Adim.';
    end
    
    fseek(fid,ElecOff+128*vOrder(iChan)+44,'bof');
    header.elec(iChan).FsCoeff = fread(fid,1,'ushort');
    fseek(fid,ElecOff+128*vOrder(iChan)+90,'bof');
    header.elec(iChan).XPos = fread(fid,1,'float');
    header.elec(iChan).YPos = fread(fid,1,'float');
    header.elec(iChan).ZPos = fread(fid,1,'float');
    fseek(fid,ElecOff+128*vOrder(iChan)+102,'bof');
    header.elec(iChan).Type = fread(fid,1,'ushort');
end
header.elec = header.elec;

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
  % determine the header.election of data to read
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

%   % output the data
  for iElec = 1 : header.Num_Chan
         data(iElec,:) = ((data(iElec,:)-header.elec(iElec).LogicGnd)/(header.elec(iElec).LogicMax-header.elec(iElec).LogicMin+1)) ...
        *(header.elec(iElec).PhysMax-header.elec(iElec).PhysMin);
  end;
  output = data;
  % FIXME why is this value of -32768 subtracted?
  % FIXME some sort of calibration should be applied to get it into microvolt
  % FIXED by Romain Ligneul / 10-2-2015
end



fclose(fid);

