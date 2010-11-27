function ds=readCTFds(datasetname)

% ************************************************************************
%
%   This program is provided to users of CTF MEG systems as a courtesy only.
%   It's operation has not been verified for clinical use.
%   Please do not redistribute it without permission from CTF Systems Inc.
%
% ************************************************************************

%  readCTFds opens a CTF data set and creates a MATLAB structure with the information
%  in the .res4, .infods, .newds, .acq, .hc, .hist, .EEG, bad.segments, BadChannels,
%  ClassFile.cls andMarkerFile.mrk files.  It confirms the existence of a .meg4 file
%  with the correct size.

%  See document CTF MEG File Formats, PN900-0088 for a description of the formats of
%  dataset files.

%  Author : Harold Wilson

% *****************  Revisions and bug fixes  **************************************
%  Version 1.3  4 October 2007
%  1. readCTFds v1.2 failed when run in MATLAB 7.2, but ran correctly in
%    MATLAB 6.5 and 7.3.  The failure occurred when calls to fscanf with a
%    final '\n' in the format string were followed by fgetl.  In MATLAB 6.5 and 7.3,
%    this returns the next line of a text file, but in MATLAB 7.2, this
%    returns an empty charatcer string.
%    Changes were made to subprograms
%           - readHc,
%           - readClassFile
%           - readMarkerFile.
%  2. In v1.1 and 1.2, if any of the head coil positions exceeded 100 cm, readHc
%     reported an error and returned hc=struct([]).  In v1.3, is reports the error,
%     and returns the erroneaous values. 


%  Version 1.2:  24 April 2007
%          - readHc modified to read fMEG .hc files.
%          - readCPersist modified to handle extra bytes that appear in some .acq files.

%  Version 1.1:  13 April 2007
%  readCTFds is modified to handle cases where the data files exceed a total of 2^31-8
%  bytes.  In these cases files .1_meg4,.2_meg4,... appear in the dataset.
% ***********************************************************************************

%  Inputs : datasetname : Complete name of the data set directory.  Includes the complete
%                         path and the .ds extension.

%  Output: ds : A structure that gives data set parameters.

%  Function calls included in this listing:
%     - readRes4
%     - readHc
%     - readEEG
%     - readMarkerFile
%     - readClassFile
%     - readBadSegments
%     - readVirtualChannels
%  External functions: - readCPersist    Reads infods and acq files.

%    - The data directory is datasetname.
%    - path is the complete path to the directory including the last file delimiter.
%    - baseName is the directory name, less the last file extension.
%    - datasetname=[path,baseName,'.ds'];

persistent printWarning multipleMeg4Files

if nargin==0 & nargout==0  %  Print a version number
  fprintf(['\treadCTFds: Version 1.3   4 October 2007   ',...
    'Reads v4.1 and v4.2 CTF data sets including fMEG .hc files.\n',...
    '\tCall: ds=readCTFds(datasetname)\n',...
    '\t\tdatasetname = Name of the dataset including the path.\n',...
    '\t\tds = Structure containing all dataset information except for data in ',...
    'the .meg4 file.\n\n']);
  return
end

MAX_COILS=8;
MAX_BALANCING=50; % max. no. of balancing coefficients
SENSOR_LABEL=31;  % No. of characters in sensor labels in the balancing
%                   coefficient tables
len_sensor_name=32;

% length of sensor coefficient records in bytes.  See listing for MEGDefs.h
%senres_lenrec=len_sensor_name+8+2+MAX_BALANCING*SENSOR_LABEL+8*MAX_BALANCING;

%  Allowed 8-byte headers for res4 and meg4 files.
res4Headers=strvcat(['MEG41RS',char(0)],['MEG42RS',char(0)]);
meg4Headers=strvcat(['MEG41CP',char(0)],['MEG42CP',char(0)]);

delim=filesep;
ds=struct([]);

%  Check that the data set exists.
if exist('datasetname')~=1
  fprintf('\nreadCTFds: No input datasetname specified.\n\n');
  return
elseif ~ischar(datasetname) | isempty(datasetname) | size(datasetname,1)~=1
  fprintf('\nreadCTFds: Dataset name is not a string, or is empty, or is an array.\n\n');
  whos datasetname
  return
else
  %  Separate datasetname into a path and the baseName and add extension .ds.
  datasetname=deblank(datasetname);
  ksep=max([0 findstr(datasetname,delim)]);
  baseName=datasetname((ksep+1):length(datasetname));
  path=datasetname(1:ksep); % String path is terminated by the file delimiter (or path=[]).
  %  Remove the last .ds from baseName.
  kdot=max(findstr(baseName,'.ds'));
  if kdot==(length(baseName)-2)
    baseName=baseName(1:(max(kdot)-1));
  else
    datasetname=[datasetname,'.ds'];
  end
  clear ksep kdot;
  if exist(datasetname)~=7
    fprintf('\nreadCTFds: Cannot find dataset %s.\n\n',datasetname);
    return
  end
end

%  Check that the res4 and meg4 files exist.
res4File=[datasetname,delim,baseName,'.res4'];
meg4File=[datasetname,delim,baseName,'.meg4'];

if exist(res4File)~=2 | exist(meg4File)~=2
  fprintf('readCTFds: In directory %s, cannot find .res4 and/or .meg4 files.\n',...
    datasetname);
  return
end

%  Check the headers on .meg4, .1_meg4, ...
qFile=0;
meg4Ext='.meg4';
meg4Size=[];
while 1
  if qFile>0;
    meg4Ext=['.',int2str(qFile),'_meg4'];
    meg4File=[datasetname,delim,baseName,meg4Ext];
  end
  if qFile>1 & isempty(multipleMeg4Files)
    fprintf('readCTFds: This dataset has multiple meg4 files.\n');
    multipleMeg4Files=1;
  end
  fid=fopen(meg4File,'r','ieee-be');
  if fid<=0;break;end
  D=dir([datasetname,delim,baseName,meg4Ext]);
  meg4Size=[meg4Size D.bytes];
  meg4Header=char(fread(fid,8,'uint8')');
  fclose(fid);
  if isempty(strmatch(meg4Header,meg4Headers,'exact'))
    fprintf('\nreadCTFds: %s file header=%s   Valid header options:',meg4Ext,meg4Header);
    for k=1:size(meg4Headers,1);fprintf('  %s',meg4Headers(k,:));end
    fprintf('\n\n');
    ds=struct([]);
    return
  end
  qFile=qFile+1;
end
Nmeg4=length(meg4Size);
clear D fid qFile;
% Add baseName and path to structure ds.
ds=struct('baseName',baseName,'path',path);

%  Read the res4 file
ds.res4=readRes4(res4File,res4Headers,...
  MAX_COILS,MAX_BALANCING,SENSOR_LABEL,len_sensor_name);
if isempty(ds.res4);ds=struct([]);return;end

dataBytes=4*ds.res4.no_trials*ds.res4.no_channels*ds.res4.no_samples;

%  Assemble ds.meg4.
ds.meg4.fileSize=sum(meg4Size);
ds.meg4.header=meg4Header;
clear meg4Header meg4Size;

%  Does the size of the .meg4 file match the size specified by the .res4 file?
if ds.meg4.fileSize~=8*Nmeg4+dataBytes
  fprintf(['\nreadCTFds: Data set error : size of meg4 file(s)\n\t\t',...
    '%10d bytes (from dir command)\n'],ds.meg4.fileSize);
  fprintf('\t\t%10d bytes (from res4 file)\n\n',8*Nmeg4+dataBytes);
  return
end

if isempty(printWarning)
  fprintf(['\nreadCTFds: You are reading CTF data for use with a software-application tool\n',...
    '\tthat is not manufactured by VSM MedTech Ltd. and has not received marketing\n',...
    '\tclearance for clinical applications.  If CTF MEG data are processed by this tool,\n',...
    '\tthey should not be later employed for clinical and/or diagnostic purposes.\n\n']);
  printWarning=1;
end

% .infods file
if exist([datasetname,delim,baseName,'.infods'])==2
  ds.infods=readCPersist([datasetname,delim,baseName,'.infods']);
end

% .newds file
if exist([datasetname,delim,baseName,'.newds'])==2
  fid=fopen([datasetname,delim,baseName,'.newds'],'r','ieee-be');
  ds.newds=char(fread(fid,'uint8'))';
  fclose(fid);
  clear fid;
end

%  .acq file
if exist([datasetname,delim,baseName,'.acq'])==2
  ds.acq=readCPersist([datasetname,delim,baseName,'.acq']);
end

%  .hist file
if exist([datasetname,delim,baseName,'.hist'])==2
  fid=fopen([datasetname,delim,baseName,'.hist'],'r','ieee-be');
  ds.hist=char(fread(fid,'uint8'))';
  fclose(fid);
end

%  .hc file
if exist([datasetname,delim,baseName,'.hc'])==2
  ds.hc=readHc([datasetname,delim,baseName,'.hc']);
end

%  .eeg file
if exist([datasetname,delim,baseName,'.eeg'])==2
  ds.eeg=readEEG([datasetname,delim,baseName,'.eeg']);
  if isempty(ds.eeg);ds=rmfield(ds,'eeg');end
end

%  marker file
if exist([datasetname,delim,'MarkerFile.mrk'])==2
  ds.mrk=readMarkerFile([datasetname,delim,'MarkerFile.mrk']);
end

%  ClassFile
if exist([datasetname,delim,'ClassFile.cls'])==2
  ds.TrialClass=readClassFile([datasetname,delim,'ClassFile.cls']);
end

%  bad.segments
if exist([datasetname,delim,'bad.segments'])==2
  ds.badSegments=readBadSegments([datasetname,delim,'bad.segments']);
end

%  BadChannels
if exist([datasetname,delim,'BadChannels'])==2
  fid=fopen([datasetname,delim,'BadChannels'],'r','ieee-be');
  ds.BadChannels=char([]);
  while 1
    strng=fscanf(fid,'%s\n',1);
    if isempty(strng);break;end
    ds.BadChannels=strvcat(ds.BadChannels,strng);
  end
  fclose(fid);
  clear fid;
end

% VirtualChannels  Assume this is the name of the file.  Modify to accept *.vc?
if exist([datasetname,delim,'VirtualChannels'])==2
  ds.Virtual=readVirtualChannels([datasetname,delim,'VirtualChannels']);
end

% processing.cfg
if exist([datasetname,delim,'processing.cfg'])==2
  fid=fopen([datasetname,delim,'processing.cfg']);
  ds.processing=char(fread(fid,'uint8'))';
  fclose(fid);
  clear fid;
end

return
%%%%%%%%%%%%%%%%%     End of readCTFds     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function readRes4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res4=readRes4(res4File,res4Headers,...
  MAX_COILS,MAX_BALANCING,SENSOR_LABEL,len_sensor_name);

% Open the res4 file
fid=fopen(res4File,'r','ieee-be');

%  Check header.
res4.header=char(fread(fid,8,'uint8')');
if isempty(strmatch(res4.header,res4Headers,'exact'))
  fprintf(['\nreadCTFds (readRes4):res4 file header = %s.  ',...
    'Valid header options:'],res4.header);
  for k=1:size(res4Headers,1);fprintf(' %s',res4Headers(k,:));end;fprintf('\n\n');
  res4=struct([]);
  fclose(fid);
  return
end

% Remove trailing blanks, but include a final null (char(0)).
res4.appName=[deblank(char(fread(fid,256,'uint8')')) char(0)];
res4.dataOrigin=[deblank(char(fread(fid,256,'uint8')')) char(0)];
res4.dataDescription=[deblank(char(fread(fid,256,'uint8')')) char(0)];
res4.no_trials_avgd=fread(fid,1,'int16');
res4.data_time=[deblank(char(fread(fid,255,'uint8')')) char(0)];
res4.data_date=[deblank(char(fread(fid,255,'uint8')')) char(0)];

% new_general_setup_rec_ext part of meg41GeneralResRec
res4.no_samples=fread(fid,1,'int32');
temp=fread(fid,2,'int16');
res4.no_channels=temp(1);
res4.sample_rate=fread(fid,1,'float64');
res4.epoch_time=fread(fid,1,'float64');
temp=fread(fid,2,'int16');
res4.no_trials=temp(1);
res4.preTrigPts=fread(fid,1,'int32');
res4.no_trials_done=fread(fid,1,'int16');
res4.no_trials_display=fread(fid,1,'int16');
res4.save_trials=fread(fid,1,'int32');

% meg41TriggerData part of new_general_setup_rec_ext   10 bytes total
res4.primaryTrigger=fread(fid,1,'int32');
res4.triggerPolarityMask=fread(fid,1,'int32');
temp=fread(fid,1,'int16');            %Skip 2 bytes of padding
% end of meg41TriggerData part of new_general_setup_rec_ext

temp=fread(fid,3,'int16');            %Skip 2 bytes of padding
res4.trigger_mode=temp(2);
res4.accept_reject_Flag=fread(fid,1,'int32');

temp=fread(fid,2,'int16');            %Skip 2 bytes of padding
res4.run_time_display=temp(1);
res4.zero_Head_Flag=fread(fid,1,'int32');
res4.artifact_mode=fread(fid,1,'int32');
%  end of new_general_setup_rec_ext part of meg41GeneralResRec

% meg4FileSetup part of meg41GeneralResRec
% Remove trailing blanks, but include a final null (char(0))
res4.nf_run_name=[deblank(char(fread(fid,32,'uint8')')) char(0)];
res4.nf_run_title=[deblank(char(fread(fid,256,'uint8')')) char(0)];
res4.nf_instruments=[deblank(char(fread(fid,32,'uint8')')) char(0)];
res4.nf_collect_descriptor=[deblank(char(fread(fid,32,'uint8')')) char(0)];
res4.nf_subject_id=[deblank(char(fread(fid,32,'uint8')')) char(0)];
res4.nf_operator=[deblank(char(fread(fid,32,'uint8')')) char(0)];
res4.nf_sensorFileName=[deblank(char(fread(fid,56,'uint8')')) char(0)];
temp=fread(fid,3,'int32');
res4.rdlen=temp(2);

res4.run_description=[deblank(char(fread(fid,res4.rdlen,'uint8')')) char(0)];
%  end of meg4FileSetup part of meg41GeneralResRec

%  Filter descriptions.  Set field res4.filters=[] if no filters are
%  defined.
res4.num_filters=fread(fid,1,'int16');
if res4.num_filters==0
  res4.filters=[];
else
  for kfilt=1:res4.num_filters
    res4.filters(kfilt).freq=fread(fid,1,'float64');
    res4.filters(kfilt).fClass=fread(fid,1,'int32');
    res4.filters(kfilt).fType=fread(fid,1,'int32');
    res4.filters(kfilt).numParam=fread(fid,1,'int16');
    for kparm=1:res4.filters(kfilt).numParam
      res4.filters(kfilt).Param(kparm)=fread(fid,1,'float64');
    end
  end
  clear kfilt kparm;
end
%  Read channel names.   Must have size(res4.chanNames)=[channels 32]
%  Clean up the channel names.  The MEG system software leaves junk
%  bytes in the channel name following the first zero.

res4.chanNames=char(zeros(res4.no_channels,32));
for kchan=1:res4.no_channels
  xname=strtok(char(fread(fid,32,'uint8')'),char(0));
  res4.chanNames(kchan,1:length(xname))=xname;
end
clear kchan xname;

%  Read sensor resource table.  Floating point values are 'double'
%  but the code could be changed to convert to 'single' to save memory.
for kchan=1:res4.no_channels
  res4.senres(kchan).sensorTypeIndex=fread(fid,1,'int16');
  res4.senres(kchan).originalRunNum=fread(fid,1,'int16');
  res4.senres(kchan).coilShape=fread(fid,1,'int32');
  res4.senres(kchan).properGain=fread(fid,1,'double');
  res4.senres(kchan).qGain=fread(fid,1,'double');
  res4.senres(kchan).ioGain=fread(fid,1,'double');
  res4.senres(kchan).ioOffset=fread(fid,1,'double');
  res4.senres(kchan).numCoils=fread(fid,1,'int16');
  numCoils=res4.senres(kchan).numCoils;
  temp=fread(fid,3,'int16');
  res4.senres(kchan).grad_order_no=temp(1);

  %  Special code to take care of situations where someone wants to label bad channels
  %  by setting their gain to zero.
  invgain=(res4.senres(kchan).ioGain*...
    res4.senres(kchan).properGain*res4.senres(kchan).qGain);
  if abs(invgain)>1e-50
    res4.senres(kchan).gain=1/invgain;
  else
    res4.senres(kchan).gain=sign(invgain)*1e50;
  end
  if res4.senres(kchan).sensorTypeIndex>=0 & res4.senres(kchan).sensorTypeIndex<=7
    % Nominal gain (fT/integer step)
    res4.senres(kchan).gain=1e15*res4.senres(kchan).gain;
  end

  % Data that was included in res4.senres(kchan).coilTbl in earlier versions of readCTFds

  res4.senres(kchan).pos0=zeros(3,res4.senres(kchan).numCoils);
  res4.senres(kchan).ori0=zeros(3,res4.senres(kchan).numCoils);
  res4.senres(kchan).area=zeros(1,res4.senres(kchan).numCoils);
  res4.senres(kchan).numturns=zeros(1,res4.senres(kchan).numCoils);
  for qx=1:numCoils
    buff=fread(fid,8,'double');
    res4.senres(kchan).pos0(:,qx)=buff(1:3);
    res4.senres(kchan).ori0(:,qx)=buff(5:7);
    temp=fread(fid,4,'int16');
    res4.senres(kchan).numturns(qx)=temp(1);
    res4.senres(kchan).area(qx)=fread(fid,1,'double');
  end
  if numCoils<MAX_COILS   % Skip the rest of the coilTbl
    buff=fread(fid,10*(MAX_COILS-numCoils),'double');
  end

  % Data that was included in res4.senres(kchan).HdcoilTbl in earlier versions of readCTFds

  res4.senres(kchan).pos=zeros(3,res4.senres(kchan).numCoils);
  res4.senres(kchan).ori=zeros(3,res4.senres(kchan).numCoils);
  for qx=1:numCoils
    buff=fread(fid,8,'double');
    res4.senres(kchan).pos(:,qx)=buff(1:3);
    res4.senres(kchan).ori(:,qx)=buff(5:7);
    temp=fread(fid,2,'double');  % Don't bother with area and numturns.  Already read.
  end
  if numCoils<MAX_COILS      % Skip the rest of the HdcoilTbl
    buff=fread(fid,10*(MAX_COILS-numCoils),'double');
  end
end
clear kchan buff numCoils temp qx;
%  End reading sensor resource table

res4.numcoef=fread(fid,1,'int16');  % Number of coefficient records

%  Create structure array res4.scrr which holds the balancing tables.
for nx=1:res4.numcoef
  res4.scrr(nx).sensorName=uint8(fread(fid,len_sensor_name,'uint8'))';
  buff=fread(fid,8,'uint8');
  res4.scrr(nx).coefType=uint8(buff(1:4))';  % discard bytes 5-8
  res4.scrr(nx).numcoefs=double(fread(fid,1,'int16'));
  numcoef=res4.scrr(nx).numcoefs;
  buff=fread(fid,SENSOR_LABEL*MAX_BALANCING,'uint8');
  buff=reshape(buff,SENSOR_LABEL,MAX_BALANCING);
  res4.scrr(nx).sensor=[uint8(buff(:,1:numcoef)) ...
    uint8(zeros(SENSOR_LABEL,MAX_BALANCING-numcoef))];
  buff=fread(fid,MAX_BALANCING,'double');
  res4.scrr(nx).coefs=[buff(1:numcoef);zeros(MAX_BALANCING-numcoef,1)];
end
fclose(fid); % Close the res4 file
return

%%%%%%%%%%  End of readRes4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function readHc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hc=readHc(hcFile);

%  Reads the .hc file and returns a structure with head-coil posiiotns.
%  The hc file format is given in document CTF MEG File Formats (PN900-0088).

%  .hc file format: Coils and positions are specified by 4 lines of text:
%  Line 1: A B coil position relative to C (cm):
%    "  2: [tag]x = F
%    "  3: [tag]y = F
%    "  4: [tag]z = F

%    A = 'standard' or 'stadard' or 'measured'.  'stadard' is a typographical error that
%        appeared in v4 of the CTF Acq.  The first set of coil positions is 'standard' or
%        'stanard'.  The next two sets are 'measured'.
%    B = coil label.  For MEG systems, the first three coils must have the names 'nasion',
%        'left ear', 'right ear'.  They may be followed by other coil names.
%         For fMEG systems, there is no fixed standard for coil names.b  readHc will
%         accept any coil names.  Typically the coils are 'back','left hip','right hip',
%         'sternum', 'belly'.  In some UAMS data the names are 'nasion',...

%        The string 'coil position relative to' appears in all valid .hc files.  It marks
%        the end of the coil label and the start ot the coordinate-system.

%    C = coordinate system = 'dewar' and 'head'(MEG) or 'abdomen' (fMEG).  No variations
%        are permitted.  The first two sets of coil positions are in dewar coordinates.
%        The last set is in 'head' or 'abdomen' coordinates.

%    Units : For .hc to be a valid coil position file, the coordinates must be in cm.
%            The program searches for the string '(cm):'.

%    F = coordinate (numeric)

%   There are exactly 3 sets of coil positions in the file appearing in the order
%          (1) A='standard', C='dewar'
%          (2) A='measured', C='dewar'
%          (3) A='measured', C='head' or 'abdomen'
%   The coils must appear in the same order in each set.

%  Input : Name of hc file (including complete path and .hc extension)

%  Output : Structure hc containing standard, dewar and CTF coordinates of the
%           nasion,left,right coils.
%             hc.name= names of the coils  size(hc.name)=[nCoil x]
%      MEG: hc.standard, hc.dewar, hc.head : size=[3 nCoil] : coil positions in cm
%     fMEG: hc.standard, hc.dewar, hc.abdomen : size=[3 nCoil] : coil positions in cm

if nargin==0 & nargout==0
  fprintf(['readHc:  Version 1.3   4 Oct. 2007\n',...
    '\thc=readHc(hcFile) reads head-coil file hcFile and returns the head coil',...
    '\t   positions in structure hc\n\n',...
    '\tThe file format is defined in document CTF MEG File Formats, PN900-0088.\n\n']);
  return
end

basicString='coil position relative to';
%  In MEG systems, the first three coil names are fixed.
MEGCoilLabel=strvcat('nasion','left ear','right ear');
C3options=strvcat('head','abdomen');
unitString='(cm):';
maxRCoil=200;  % Report an error if any coil position is > maxRCoil cm.
printWarning=1;
hc=struct([]);  % Return in event of an error.

if exist(hcFile)~=2
  return
end

prec=1e-6;  % Round positions.
%  Decide if this is a MEG or an fMEG system.  Find the first line that contains the
%  string 'coil position relative to'.
%        MEG system: This line contains 'nasion'.
%       fMEG system: This line does not contains 'nasion'.

fhc=fopen(hcFile,'rt','ieee-be');

%  Decide system type and generate arrays for strings A and C.
%  Search for first line.  Allow for the posibility that extra text lines migt be added at
%  the start of the file.
strngA=char([]);

while 1
  textline=fgetl(fhc);
  if isequal(textline,-1)
    break
  elseif strfind(textline,'coil position relative to');
    stringA=strtok(textline);
    % Fix typographic error in some releases of Acq, not in all coils.
    if strcmp(stringA,'stadard')
      stringA = 'standard';
    end
    if ~strcmp(stringA,'standard')
      break
    end
    stringA=strvcat(stringA,'measured','measured');
    break;
  end
end

if isequal(textline,-1) | size(stringA,1)<3
  fprintf('readHc: File %s does not have the head position file format.\n',hcFile);
  fclose(fhc);
  return
end

%  Set stringC(3,:)='head' or 'abdomen' later.
stringC=strvcat('dewar','dewar');

%  textline is the first useful line of the head coil file.

coilLabel=char([]);   % Names given to the coils
rhc=zeros(3,0);       % Coil positions (cm)
coil=0;
nCoil=0;
lastPositionSet=1;

while ~isequal(textline,-1)
  %  Parse the line
  [A,R]=strtok(textline);
  % Fix typographic error in some releases of Acq, not in all coils.
  if strcmp(A,'stadard')
    A = 'standard';
  end
  kpos=strfind(R,basicString);
  coilName=deblank(R(2:kpos-1));
  [C,R]=strtok(R(kpos+length(basicString):length(R)));
  if strmatch(C,C3options)
    stringC=strvcat(stringC,C);
  end
  [unit,R]=strtok(R);
  positionSet=intersect(strmatch(A,stringA),strmatch(C,stringC));
  if isempty(positionSet) | ~strcmp(unit,unitString) | ~isempty(R)
    break;
  end
  if positionSet==lastPositionSet
    coil=coil+1;
  else
    coil=1;
  end
  if positionSet==1
    %  Assemble list of coil names
    coilLabel=strvcat(coilLabel,coilName);
    nCoil=coil;
  elseif ~strcmp(coilName,deblank(coilLabel(coil,:)));
    break;
  end
  %  The line describing the coil and coordinate frame is OK.
  %  Get the coordinates by reading 3 lines.
  buff=fscanf(fhc,'%s%s%f',9);
  if ~isequal(size(buff),[9 1])
    break
  end
  rhc=[rhc prec*round(buff(3:3:9)/prec)];
  fgetl(fhc);  % Skip to the start of the next line.
  textline=fgetl(fhc);
  lastPositionSet=positionSet;
end
fclose(fhc);

clear lastPositionSet coil textline fhc;
clear stringA A C R kpos coilName unit buff;

if size(rhc,2)~=positionSet*nCoil
  fprintf('readHc: File %s does not have %d complete sets of %d coils.\n',...
    hcFile,positionSet,nCoil);
  return
end

if max(max(abs(rhc)))>maxRCoil & printWarning
  fprintf('readHc: Head Coil file %s\n        max(coil position)>%d.\n',...
    hcFile,round(maxRCoil));
end

%  Assemble structure hc.
hc=struct('names',coilLabel);
coilCoords=deblank(strvcat('standard','dewar',stringC(3,:)));
for q=1:positionSet
  hc=setfield(hc,deblank(coilCoords(q,:)),rhc(:,nCoil*(q-1)+[1:nCoil]));
end
return
%%%%%%%%%  End of readHc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function readEEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EEG=readEEG(EEGfile);

%  Reads the EEG file of a dataset and stores the infoemation in strucure array EEG.
%  EEG(k).chanName = channel name in the dataset (EEGmmm where mmm=channel number)
%  EEG(k).name = channel name given by the user (e.g. Fp4)
%  EEG(k).pos = electrode position in cm in CTF head coordinates

if exist(EEGfile)~=2
  EEG=struct([]);
  return
end

fid=fopen(EEGfile,'r');
EEG=struct('chanNum',[],'name',char([]),'pos',[]);
nEEG=0;
while 1
  chanNum=fscanf(fid,'%d',1);
  name=fscanf(fid,'%s',1);
  pos=fscanf(fid,'%f',3);
  if isempty(pos);break;end
  nEEG=nEEG+1;
  EEG(nEEG)=struct('chanNum',chanNum,'name',name,'pos',pos);
end
fclose(fid);
return
%%%%%%%%%  End of readEEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function readClassFile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TrialClass=readClassFile(ClassFileName);

%   Reads a CTF ClassFile and stores information in a structure.
%   The class file allows a user to store a list of classsified trials in a data set.
%   The ClassFile format is defined in document CTF MEG File Formats, PN900-0088.
%   This format is rigid and readClassFile assumes that the ClassFile has the format
%   current in October 2006.

%   Inputs :
%      ClassFileName : marker file including the full path and extension .mrk.
%      trialList : List of trials to read.  Trial numbering : 1,2,...
%                  If omitted or empty, read all markers.

%  Output : Structure array marker.  Output trial numbering starts at 1.
%           See CTF MEG File Formats, (PN900-0088) for the meaning of the structure
%           fields.  A trial mat=y start before the t=0 point, so it is possible to have
%           markers with time<0 (see ds.res4.preTrigPts).

TrialClass=struct([]);

if exist(ClassFileName)~=2
  return     % File doesn't exist.
end

fid=fopen(ClassFileName,'r','ieee-be');

for k=1:5;fgetl(fid);end  % Skip 5 lines (including path info)
nClass=sscanf(fgetl(fid),'%d',1); %Read value and skip to the start of the next non-blank line.
if nClass<=0
  fprintf('readClassFile: File %s has %d classes.\n',nClass);
  return
end

TrialClass=struct('ClassGroupId',[],'Name',char([]),...
  'Comment',char([]),'Color',char([]),'Editable',char([]),'ClassId',[],'trial',[]);

for k=1:nClass
  %  Find the start of the next class identification
  %  There is no need to check for end of file because the loop ends before an attempt
  %  is made to read class nClass+1.
  while ~strcmp('CLASSGROUPID:',fgetl(fid));end
  ClassGroupId=sscanf(fgetl(fid),'%d',1);
  fgetl(fid);
  Name=deblank(fgetl(fid));
  fgetl(fid);
  Comment=deblank(fgetl(fid));
  fgetl(fid);
  Color=deblank(fgetl(fid));
  fgetl(fid);
  Editable=deblank(fgetl(fid));
  fgetl(fid);
  ClassId=sscanf(fgetl(fid),'%d',1);
  fgetl(fid);
  No_of_Trials=sscanf(fgetl(fid),'%d',1);
  fgetl(fid);fgetl(fid);
  if No_of_Trials>0
    trial=reshape(fscanf(fid,'%d',No_of_Trials),1,No_of_Trials);
  else
    trial=[];
  end
  %  Adjust trial numbering so it starts at 1.
  TrialClass(k)=struct('ClassGroupId',ClassGroupId,'Name',Name,...
    'Comment',Comment,'Color',Color,'Editable',Editable,'ClassId',ClassId,...
    'trial',trial+1);
end
fclose(fid);
return
%%%%%%%%%  End of readClassFile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function readBadSegments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function badSegments=readBadSegments(badSegmentsFile);

%  Reads the bad.segements file of a CTF data set and stores the information in structure
%  bad_segments.
%    badSegments.trial = List of trial numbers
%    badSegments.StartTime = List of bad segment start times (relative to trial).
%    badSegments.EndTime   = List of bad segment end times.

%  Reads the file one line at a time.  Each line has three numbers.
%  1st: Trial number=integer.  Trial numbering for bad.segments starts at 1.
%  2nd,3rd : Start, End times (s).

if exist(badSegmentsFile)~=2
  badSegments=struct([]);
  return
end

fid=fopen(badSegmentsFile,'r','ieee-be');
badSegments=struct('chanName',char([]),'name',char([]),'pos',[]);
nLine=0;
data=zeros(3,0);
while 1
  txt=fgetl(fid);
  if isequal(txt,-1);break;end
  nLine=nLine+1;
  [buff,count]=sscanf(txt,'%f');
  %  Are the data good?
  badline=(count~=3);
  if ~badline;badline=(abs((buff(1)-round(buff(1)))>0.0001) | buff(3)<buff(2));end
  if badline
    fprintf('readCTFds (readBadSegments): Format error in file bad.segments.\n');
    fprintf('\tLine %d has %d numbers: ',nLine,count);fprintf('\t %g',buff);fprintf('\n');
    fprintf('\tMust have 3 numbers on each line: #1=integer, #2,#3=float, #3>=#2.\n');
    badSegments=struct([]);
    return
  end
  data=[data buff];
end
fclose(fid);
badSegments=struct('trial',data(1,:),'StartTime',data(2,:),'EndTime',data(3,:));

return
%%%%%%%%%  End of readBadSegments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function readVirtualChannels  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Virtual=readVirtualChannels(VirtualChannelsName);

%  Reads the Virtual channels file associated with a data set

Virtual=struct([]);
%  Check VirtualChannelsName
if exist('VirtualChannelsName')~=1;VirtualChannelsName=char([]);end
if isempty(VirtualChannelsName)
  fprintf('readVirtualChannels:  Must specify a VirtualChannels file.\n');
  return
elseif ~ischar(VirtualChannelsName)
  fprintf('readVirtualChannels:  VirtualChannelsName must be a character string.\n');
  return
elseif exist(VirtualChannelsName)~=2
  fprintf('readVirtualChannels:  Cannot find file %s\n',VirtualChannelsName);
  return
end

fid=fopen(VirtualChannelsName,'r');
count=0;
Virtual=struct('Name',char([]),'Unit',char([]),'chan',char([]),'wt',[]);

strng=textread(VirtualChannelsName,'%s','delimiter',',\t');
k=0;
while k<length(strng)
  k=k+1;
  if strmatch('VirtualChannel',strng{k});
    Refcount=0;
    chan=char([]);
    wt=[];
  elseif strmatch('Name:',strng{k});
    k=k+1;
    Name=strng{k};
    if isempty(Name);break;end
  elseif strmatch('Unit:',strng{k});
    k=k+1;
    Unit=strng{k};
  elseif strmatch('Ref:',strng{k});
    chan=strvcat(chan,strng{k+1});
    wt=[wt;str2num(strng{k+2})];
    k=k+2;
  elseif strmatch('}',strng{k});
    count=count+1;
    Virtual(count)=struct('Name',Name,'Unit',Unit,'chan',chan,'wt',wt);
    clear Name Unit chan wt Refcount;
  end
end
fclose(fid);
if isempty(Virtual(1).Name)
  Virtual=struct([]);
end
return
%%%%%%%%%  End of readVirtualChannels  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function readMarkerFile  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function marker=readMarkerFile(MarkerFileName,trialList);

% Version 1.3   4 Oct. 2007

%   Reads specified trials of a CTF MarkerFile.
%   The MarkerFile format is defined in document CTF MEG File Formats, PN900-0088.
%   This format is rigid and readMarkerFile assumes that the MarkerFile has the format
%   current in October 2006.    The document lists 4 Class Groups, but only CLASSGROUPID=0
%   (triggers) and CLASSGROUPID=3 (Manual) are actually used.  read_MArkerFile reads only
%   these two groups, but it prints a message if the other groups are encountered.

%   Trigger markers (CLASSGROUPID=0) have an additional 4 pieces of information supplied:
%    BITNUMBER(int), POLARITY(+ or 1),SOURCE(text) and THRESHOLD (float).  FOr
%    CLASSGROUPID=3, these fields are left empty.


%   Inputs :
%      MarkerFileName : marker file including the full path and extension .mrk.
%      trialList : List of trials to read.  Trial numbering : 1,2,...
%                  If omitted or empty, read all markers.

%  Output : marker: Structure array. Output trial numbering starts at 1.
%            See CTF MEG File Formats, (PN900-0088) for the meaning of the structure
%            fields.  A trial may start before the t=0 point, so it is possible to have
%            markers with time<0 (see ds.res4.preTrigPts).
%            Fields: ClassGroupId, Name,Comment,Color,Editable,ClassId,trial,time;
%                    Other fields that appear with ClassGroupId=0 are not included.


if nargin==0 & nargout==0
  fprintf(['readMarkerFile:  Version 1.3   4 Oct. 2007\n',...
    '\tReads a CTF dataset MarkerFile and returns a structure array.\n',...
    '\tmarker=readMarkerFile(FileName) reads ',...
    'a marker file and returns the trials and times.\n',...
    '\tmarker=readMarkerFile(FileName,trialList) reads a marker file and returns only\n',...
    '\t\t\tthe markers for trials listed in vector trialList.\n\n',...
    '\treadMarkerFile increments trial numbers by 1 so first trial in a dataset has trial=1.\n\n',...
    '\tThe MarkerFile format is defined in document CTF MEG File Formats, PN900-0088.\n\n']);
  return
end

marker=struct([]);

%  Check MarkerFileName
if exist('MarkerFileName')~=1;MarkerFileName=char([]);end
if exist(MarkerFileName)~=2   % No MarkerFile present
  fprintf('readMarkerFile:  Could not find MarkerFile.\n');
  MarkerFileName
  return
end

%  Check input trialList
if exist('trialList')~=1
  trialList=[];
else
  trialList=reshape(trialList,1,length(trialList));
  if any((trialList-round(trialList))>0.1) | any(round(trialList)<=0)
    fprintf('readMarkerFile: trialList must have only positive integers.\n');
    return
  end
end

fid=fopen(MarkerFileName,'r','ieee-be');

for k=1:5;fgetl(fid);end  % Skip 5 lines (including path info)
nMarker=sscanf(fgetl(fid),'%d',1); %Read value and skip to the start of the next non-blank line.
if nMarker<=0
  fprintf('readMarkerFile: File %s has %d markers.\n',nMarker);
  fclose(fid);
  return
end

marker=struct('ClassGroupId',[],'Name',char([]),...
  'Comment',char([]),'Color',char([]),'Editable',char([]),'ClassId',[],...
  'BitNumber',[],'Polarity',char([]),'Source',char([]),'Threshold',[],...
  'trial',[],'time',[]);

for k=1:nMarker
  %  Find the start of the next marker identification
  %  There is no need to check for end of file because the loop ends before an attempt
  %  is made to read marker class nClass+1.
  while ~strcmp('CLASSGROUPID:',fgetl(fid));end
  ClassGroupId=sscanf(fgetl(fid),'%d',1);
  if ~any(ClassGroupId==[0 3])
    fprintf('read_MarkerFile: Skipping a marker with CLASSGROUPID=%d\n',ClassGroupId);
    continue;
  end
  fgetl(fid);
  Name=deblank(fgetl(fid));
  fgetl(fid);
  Comment=deblank(fgetl(fid));
  fgetl(fid);
  Color=deblank(fgetl(fid));
  fgetl(fid);
  Editable=deblank(fgetl(fid));
  fgetl(fid);
  ClassId=sscanf(fgetl(fid),'%d',1);
  if ClassGroupId==0
    fgetl(fid);
    BitNumber=sscanf(fgetl(fid),'%d',1);
    fgetl(fid);
    Polarity=deblank(fgetl(fid));
    fgetl(fid);
    Source=deblank(fgetl(fid));
    fgetl(fid);
    Threshold=sscanf(fgetl(fid),'%f',1);
  else
    BitNumber=[];
    Polarity=char([]);
    Source=char([]);
    Threshold=[];
  end
  fgetl(fid);
  No_of_Samples=sscanf(fgetl(fid),'%d',1);
  fgetl(fid);fgetl(fid);
  trial=zeros(1,No_of_Samples);
  time=zeros(1,No_of_Samples);
  if No_of_Samples>0
    buff=fscanf(fid,'%d %f\n',2*No_of_Samples);
    trial=reshape(buff(1:2:2*No_of_Samples-1),1,No_of_Samples)+1; % Trial numbering starts at 1.
    time=reshape(buff(2:2:2*No_of_Samples),1,No_of_Samples);
    clear buff;
  end
  %  Keep only the specified trials.
  if ~isempty(trialList)
    index=[];
    for q=trialList;
      index=[index find(trial==q)];
    end
    trial=trial(index);
    time=time(index);
    No_of_Samples=length(index);
    clear q index;
  end
  marker(k)=struct('ClassGroupId',ClassGroupId,'Name',Name,...
    'Comment',Comment,'Color',Color,'Editable',Editable,'ClassId',ClassId,...
    'BitNumber',BitNumber,'Polarity',Polarity,'Source',Source,'Threshold',Threshold,...
    'trial',trial,'time',time);
end
fclose(fid);
%clear k ClassGroupId Name Comment Color Editable ClassID No_of_Samples trial time;
return
%%%%%%%%%%%%  End of readMarkerFile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
