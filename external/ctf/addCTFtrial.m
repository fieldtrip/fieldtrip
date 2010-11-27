function cntrl=addCTFtrial(datasetname,data,unit,mrk,cls,badSegments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                %
%   This program creates datasets that can be analyzed by CTF software.          %
%                                                                                %
%   Datasets created by this program MUST NOT BE USED FOR CLINICAL APPLICATIONS. %
%                                                                                %
%   Please do not redistribute it without permission from VSM MedTech Ltd.       %
%                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  addCTFtrial.m   Version 1.1  Adds trial(s) to a a CTF-format data set previously
%         created from MATLAB by writeCTFds.m.  The purpose of addCTFtrial is to allow
%         users to write CTF data sets that are too large for the computer's memory.

%  Date : 18 April 2007
%  Author : Harold Wilson

%  Operation :
%    1.  Use writeCTFds to create a dataset (name = DATASET) consisting of one or 
%        more trials.   size(data)= [SAMPLES, CHANNELS, TRIALS].
%    2.  Later, a trial is to be added to DATASET.  This program adds the new information
%        to the .res4, .mrk, .cls and badSegements files and appends the new data to the
%        existing .meg4 file.

%  Inputs : datasetname : Dataset including path.  Extension '.ds' is optional.
%           data : MEG data array.  size(data)=[no_samples no_channels no_trials_added]
%                  Array data may be single or double.
%           unit :  Determines the unit of the SQUID and EEG signals:
%                   If unit is missing or unit==[], the unit is set to 'fT'.
%                   'ft' or 'fT' : Convert to fT (MEG), uV (EEG)
%                   't'  or 'T'  : Convert to  T (MEG),  V (EEG)
%                   'phi0' : Convert to phi0 (MEG), uV (EEG)
%                   'int': Read plain integers from *.meg4-file
%           mrk,cls,badSegments : Structures describing the markers, trial classes and bad 
%                   segments of the new trials.  Trial numbering refers to the new trials
%                   being added to the dataset, starting at trial=1.  Markers are updated
%                   only when there is already as .mrk file in the dataset.

%  Outputs : - Dataset with trial(s) added.
%            - cntrl : 0 - failure.  The dataset is not modified.
%                      1 - success. The new data is added to the dataset

%  Function calls
%   Included in this listing:
%      - combineMarkers:  Combines marker specifications in structures ds.mrk and mrk.
%      - combineClasses:  Combines trial class specifications in structures 
%                          ds.TrialClass and cls.
%      - writeBadSegments:  Creates the new bad.segments file.
%      - writeClassFile:    Creates the new ClassFile.cls file.
%      - updateCreatorSoftware:  Adds non-clinical-use and creation software messages to
%                 infods, res4, newds and hist fields of ds.
%      - updateHLC: Adds head-coil movement information from the trials being added to an
%                   existing ds.infods.
%      - getArrayField : Extracts one field of a structure array to make it easier to
%                        manipulate.

%  Not included in this listing:
%      - writeRes4: Writes the new .res4 file.
%      - writeMarkerFile:  Creates the new MarkerFile.mrk file.
%      - writeCPersist:    Creates the new .acq and .infods files.

%  Output files are opened with 'w' permission and 'ieee-be' machine format in order
%  to be compatible with the Linux acquisition and analysis software.  Do not open files
%  with 'wt' permission because this will add an extra byte (char(13)) at the end of each
%  line of text.

%persistent printWarning  % No clinical-use warning with addCTFtrial.
%printWarning=1;
delim=filesep;

if nargin==0 & nargout==0  %  Print a version number
  fprintf(['\taddCTFtrial: Version 1.1   18 April 2007   ',...
      'Adds trials to v4.1 and v4.2 CTF data sets.\n',...
      '\tCall: cntrl=addCTFtrial(datasetname,data,unit,mrk,cls,badSegments);\n',...
      '\t\tdatasetname = Name of existing dataset including the path.\n',...
      '\t\tdata = Data that will be written to the new dataset .meg4 file.\n',...
      '\t\tunit = Physical units of the data.\n',...
      '\t\tmrk,cls,badSegments = descriptions of markers, trial classes and',...
      ' bad segments of data.\n\n']);
  return
end
%  Indicate failure on an early return.
cntrl=0;

%  Allowed 8-byte headers for res4 and meg4 files.
res4_headers=strvcat(['MEG41RS',char(0)],['MEG42RS',char(0)]);
meg4_headers=strvcat(['MEG41CP',char(0)],['MEG42CP',char(0)]);

MAX_COILS=8;      % Parameter that determines the size of the ds.res4.senres structure.
maxMEG4Size=2^31;  % Maximum MEG$ file in bytes.  (Limit set by Linux software)

% addCTFds will add trials to datasets only if they were created by writeCTFds.
originalCreatorSoftware='writeCTFds';
creatorSoftware='addCTFtrial';  % Added to .infods file
meg4ChunkSize=2^20; % Write new .meg4 file in chunks of 4*meg4ChunkSize bytes.

clinical_use_message='NOT FOR CLINICAL USE';

if ~exist('clinical_use_message');
  clinical_use_message=char([]);
end

%  Check number of inputs
if nargin<2
  fprintf(['addCTFtrial: Must supply inputs datasetname,data.  ',...
      'Only %d input arguments are present.\n'],nargin);
  return
end

%  Check input argument unit.  Convert unit to lower case.
if exist('unit')~=1
  unit='ft';       % default
elseif isempty(unit)
  unit='ft';       % default
elseif ischar(unit)
  unit=lower(unit);
  if ~strcmp(unit,'int') & ~strcmp(unit,'ft') & ~strcmp(unit,'t') & ~strcmp(unit,'phi0')
    fprintf(['addCTFtrial : unit=%s    Not a valid option.  Must be ',...
        '''fT'', ''T'', ''phi0'' or ''int''\n'],unit);
    ds=-1;  % Force an error in the calling program.
    return
  end
end

% Check argument type
if ~isnumeric(data)
  fprintf('\naddCTFtrial: Input data is not numeric.\n');
  return
elseif ~ischar(unit)
  fprintf('\naddCTFtrial: Input unit is not char.\n');
  return
elseif~ischar(datasetname)
  fprintf('\naddCTFtrial: Input datasetname is not char.\n');
  return
end
if exist('mrk')==1
  if ~isstruct(mrk)
    fprintf('\naddCTFtrial: Input mrk is not a structure.\n\n');
    return
  end
else
  mrk=struct([]);
end
if exist('cls')==1
  if ~isstruct(cls)
    fprintf('\naddCTFtrial: Input cls is not a structure.\n\n');
    return
  end
else
  cls=struct([]);
end
if exist('badSegments')==1
  if ~isstruct(badSegments)
    fprintf('\naddCTFtrial: Input badSegments is not a structure.\n\n');
    return
  end
else
  badSegments=struct([]);
end

% Does the dataset exist?
if exist(datasetname)~=7
  fprintf('addCTFtrial: Dataset %s does not exist.\n',datasetname);
  return
end

%  Separate datasetname into a path and the baseName
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

%  Get the ds structure of the dataset
ds=readCTFds(datasetname);

%  Update text fields of res4,infods, newds.
ds=updateCreatorSoftware(ds,creatorSoftware,originalCreatorSoftware);
if isempty(ds);return;end

nSample=ds.res4.no_samples;
fSample=ds.res4.sample_rate;
nChan=ds.res4.no_channels;

%  Check that the new data and the existing data have the same number of channels
%  and points per trial.
if size(data,1)~=nSample | size(data,2)~=nChan | size(data,3)==0
  fprintf(['addCTFtrial: size(data)=[%d %d %d], but existing data has no_samples=%d,',...
      ' no_channels=%d\n'],[size(data,1) size(data,2) size(data,3)],...
    nSample,nChan);
  return
end

%  Update ds.res4 fields to match the size of array data.
nTrialNew=size(data,3);
nTrialOld=ds.res4.no_trials;
ds.res4.no_trials=nTrialOld+nTrialNew;
ds.res4.epoch_time=nSample*ds.res4.no_trials/fSample;

%  Before converting data to integers, save the HLC channels for motion analysis in function
%  function updateHLC.
if isempty(strmatch('HLC',ds.res4.chanNames))
  HLCdata=[];
else
  %  Make a list of head-coil channels
  coil=0;
  HLClist=[];
  while ~isempty(strmatch(['HLC00',int2str(coil+1)],ds.res4.chanNames))
    coil=coil+1;
    for k=1:3
      HLClist=[HLClist strmatch(['HLC00',int2str(coil),int2str(k)],ds.res4.chanNames)];
    end
  end
  HLCdata=reshape(double(data(:,HLClist,:)),nSample,3,coil,nTrialNew);
  clear coil k HLClist;
end

%  Convert data to integers because CTF data sets are stored as raw numbers and
%  not as physical qunatities.  The res4 file contains the calibrations for 
%  converting back to physical units.  Array data may be single precision, so 
%  convert to double before doing any adjustments to the data.
%  Already checked that unit is valid.

if strcmp(unit,'int')
  data=reshape(data,nSample,nChan*nTrialNew);
  for k=1:nChan*nTrialNew
    if strcmp(class(data),'single')
      data(:,k)=single(round(double(data(:,k))));
    else
      data(:,k)=round(double(data(:,k)));
    end
  end
  data=reshape(data,nSample,nChan,nTrialNew);
  clear k;
else
  for chan=1:nChan           %  Convert EEGs from uV to V, SQUIDs from fT to T
    SQUIDtype=any(ds.res4.senres(chan).sensorTypeIndex==[0:7]);
    EEGtype=any(ds.res4.senres(chan).sensorTypeIndex==[8 9]);
    if EEGtype & (strcmp(unit,'ft') | strtcmp(unit,'phi0'))
      alphaG=1e-6;
    elseif SQUIDtype & strcmp(unit,'ft')
      alphaG=1e-15;
    elseif SQUIDtype & strcmp(unit,'phi0')
      alphaG=1./(ds.res4.senres(chan).properGain*ds.res4.senres(chan).ioGain);
    else
      alphaG=1;
    end
    %  Convert from physical units to integers using the gains in the senres table.
    for kt=1:nTrialNew
      buff=round(double(data(:,chan,kt))*(alphaG*...
        (ds.res4.senres(chan).properGain*ds.res4.senres(chan).qGain*ds.res4.senres(chan).ioGain)));
      if strcmp(class(data),'single')
        data(:,chan,kt)=single(buff);
      else
        data(:,chan,kt)=buff;
      end
    end
  end
  clear chan alphaG SQUIDtype EEGtype buff kt;
end

%  Check the meg4 file header.

fidMeg4=fopen([path,baseName,'.ds\',baseName,'.meg4'],'r','ieee-be');
fileHeader=fread(fidMeg4,8,'char')';
fclose(fidMeg4);

if isempty(strmatch(fileHeader(1:7),meg4_headers(:,1:7),'exact'))
  ds.meg4.header=meg4_headers(1,1:7);
  fprintf('addCTFtrial: meg4 header=%s ?\n',fileHeader);
  return
end

%  Create the .res4 file in the output dataset.
ds.res4=writeRes4([path,baseName,'.ds\',baseName,'.res4'],ds.res4,MAX_COILS);

if ds.res4.numcoef<0
  fprintf('addCTFtrial: writeRes4 returned ds.res4.numcoef=%d (<0??)\n',ds.res4.numcoef);
  fprintf('             Deleting dataset %s\n',baseName);
  rmdir([path,baseName,'.ds'],'s');
  return
end

%  Add trial(s) to the existing meg4 file.  Complicated coding in case the use requests
%  writing several trials, which fill up one meg4 file and start the next one.

ptsPerTrial=nSample*nChan;
data=reshape(data,ptsPerTrial,nTrialNew);
maxTrialPerFile=floor((maxMEG4Size-8)/(4*ptsPerTrial));
maxPtsPerFile=ptsPerTrial*maxTrialPerFile;

nExt=-1;
trial=0;
pt=0;
while trial<nTrialNew
  nExt=nExt+1;
  if nExt==0
    meg4Ext='.meg4';
  else
    meg4Ext=['.',int2str(nExt),'_meg4'];
  end
  meg4File=[path,baseName,'.ds\',baseName,meg4Ext];
  D=dir(meg4File);
  if isempty(D)
    %  Create the next meg4 file and write the header.
    fidMeg4=fopen(meg4File,'w','ieee-be');
    fwrite(fidMeg4,[ds.meg4.header(1:7),char(0)],'uint8');
    fclose(fidMeg4);
    D=struct('bytes',8);
  end 
  nWrite=min(nTrialNew-trial,floor((maxMEG4Size-D.bytes)/(4*ptsPerTrial)));
  %fprintf('nExt=%d  meg4File=%s\n\t\tsize=%d bytes   nWrite=%d\n',nExt,meg4File,D.bytes,nWrite);
  if nWrite>0
    fidMeg4=fopen(meg4File,'a','ieee-be');
    endPt=pt+nWrite*ptsPerTrial;
    while pt<endPt
      pt1=min(pt+meg4ChunkSize,endPt);                 %  Convert to double in case data is
      fwrite(fidMeg4,double(data(pt+1:pt1)),'int32');  %  single.  This may cause memory
      pt=pt1;                                          %  problems so write in short pieces.
    end
    fclose(fidMeg4);
    trial=trial+nWrite;
  end
end
clear trial pt meg4Ext meg4File D fidMeg4 nWrite endPt pt1;

%  Set variable cntrl to indicate that the meg4 and res4 files have been updated.
cntrl=1;

%  Update the .meg4 part of structure ds.
ds.meg4.fileSize=8*(nExt+1);
clear data pt pt1 ndata fidMeg4;

%  Add dataset names to .hist, and write new .hist file.
if ~isfield(ds,'hist');
  ds.hist=char([]);
else
  q0=max(strfind(ds.hist,originalCreatorSoftware));
  if ~isempty(q0) & isempty(strfind(ds.hist(q0:length(ds.hist)),creatorSoftware))
    ds.hist=[ds.hist(1:q0+length(originalCreatorSoftware)-1) ', ' creatorSoftware ...
        ds.hist(q0+length(originalCreatorSoftware):length(ds.hist))];
  end
  %  Create .hist file
  histfile=[path,baseName,'.ds\',baseName,'.hist'];
  fid=fopen(histfile,'w');
  fwrite(fid,ds.hist,'char');
  fclose(fid);
end

%  Add HLC data and write new infods file.
if isfield(ds,'infods')
  ds.infods=updateHLC(ds.infods,HLCdata);
  writeCPersist([path,baseName,'.ds',delim,baseName,'.infods'],ds.infods);
end

clear HLCdata;

% New .newds file.
if isfield(ds,'newds')
  fid=fopen([path,baseName,'.ds\',baseName,'.newds'],'w');
  fwrite(fid,ds.newds,'char');
  fclose(fid);
end

% Add markers, and create new MarkerFile.mrk
if isfield(ds,'mrk')
  mrk=combineMarkers(ds.mrk,nTrialOld,mrk,nTrialNew,ds.res4);
  if ~isempty(mrk)
    writeMarkerFile([path,baseName,'.ds',delim,'MarkerFile.mrk'],mrk);
  end
elseif ~isempty(mrk)
  fprintf(['addCTFtrial: Existing dataset has no MarkerFile, but structure mrk is',...
      ' specified in call to addCTFtrial.\n',...
      '             MarkerFile.mrk not created.\n']);
end

%  bad.segments file
if isfield(ds,'badSegments') & ~isempty(badSegments)
  ds.badSegments.trial=[ds.badSegments.trial badSegments.trial+nTrialOld];
  ds.badSegments.StartTime=[ds.badSegments.StartTime badSegments.StartTime];
  ds.badSegments.EndTime=[ds.badSegments.EndTime badSegments.EndTime];
  writeBadSegments([path,baseName,'.ds',delim,'bad.segments'],ds.badSegments,...
    ds.res4.no_trials,nSample/fSample);
end

%  ClassFile.  Add trial classes, and create new ClassFile.cls
if ~isempty(cls)
  if isfield(ds,'TrialClass')
    cls=combineClasses(ds.TrialClass,nTrialOld,cls,nTrialNew);
    if ~isempty(cls)
      writeClassFile([path,baseName,'.ds',delim,'ClassFile.cls'],cls);
    end
  else
    fprintf(['addCTFtrial: Existing dataset has no ClassFile, but structure cls is',...
        ' specified in call to addCTFtrial.\n',...
        '             ClassFile.cls not created.\n']);
  end
end
return
% *************** End of function addCTFtrial ********************************************
% **************************************************************************************


% **************************************************************************************
%  ************* function make_dummy_infods*********************************************
function Tag=make_dummy_infods(exist_hc,exist_HLC,sample_rate);

%  If the user does not supply ds.infos, this code makes a dummy version.  It has all of
%  the tags that Acq created in file Richard_SEF_20060606_04.infods.
% ***********************************************************************************
% ***********************************************************************************
DATASET_HZ_UNKNOWN=round(2^31-1);  % Peculiar requirement of DataEditor as of 23 Oct. 2006

fprintf('writeCTDds (make_dummy_infods): Creating a dummy infods file with all the tags.\n');
Tag(1)=struct('name','WS1_','type',0,'data',[]);
Tag(length(Tag)+1)=struct('name','_PATIENT_INFO','type',2,'data',[]);
Tag(length(Tag)+1)=struct('name','WS1_','type',0,'data',[]);
Tag(length(Tag)+1)=struct('name','_PATIENT_UID','type',10,...
  'data','2.16.124.113000.000000.00000000000000.000000000.00000000.0011');
Tag(length(Tag)+1)=struct('name','_PATIENT_NAME_FIRST','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_PATIENT_NAME_MIDDLE','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_PATIENT_NAME_LAST','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_PATIENT_ID','type',10,'data','x');
Tag(length(Tag)+1)=struct('name','_PATIENT_BIRTHDATE','type',10,'data','19500101000000');
Tag(length(Tag)+1)=struct('name','_PATIENT_SEX','type',5,'data',2);
Tag(length(Tag)+1)=struct('name','_PATIENT_PACS_NAME','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_PATIENT_PACS_UID','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_PATIENT_INSTITUTE','type',10,'data','');
Tag(length(Tag)+1)=struct('name','EndOfParameters','type',-1,'data',[]);
Tag(length(Tag)+1)=struct('name','_PROCEDURE_INFO','type',2,'data',[]);
Tag(length(Tag)+1)=struct('name','WS1_','type',0,'data',[]);
Tag(length(Tag)+1)=struct('name','_PROCEDURE_VERSION','type',5,'data',1);
Tag(length(Tag)+1)=struct('name','_PROCEDURE_UID','type',10,...
  'data','2.16.124.113000.000000.00000000000000.000000000.00000000.0041');
Tag(length(Tag)+1)=struct('name','_PROCEDURE_ACCESSIONNUMBER','type',10,'data','0');
Tag(length(Tag)+1)=struct('name','_PROCEDURE_TITLE','type',10,'data','0');
Tag(length(Tag)+1)=struct('name','_PROCEDURE_SITE','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_PROCEDURE_STATUS','type',5,'data',1);
Tag(length(Tag)+1)=struct('name','_PROCEDURE_TYPE','type',5,'data',2);  % Research type
Tag(length(Tag)+1)=struct('name','_PROCEDURE_STARTEDDATETIME','type',10,... 
  'data','20060606164306');
Tag(length(Tag)+1)=struct('name','_PROCEDURE_CLOSEDDATETIME','type',10,...
  'data','19000100000000');
Tag(length(Tag)+1)=struct('name','_PROCEDURE_COMMENTS','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_PROCEDURE_LOCATION','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_PROCEDURE_ISINDB','type',5,'data',0);
Tag(length(Tag)+1)=struct('name','EndOfParameters','type',-1,'data',[]);
Tag(length(Tag)+1)=struct('name','_DATASET_INFO','type',2,'data',[]);
Tag(length(Tag)+1)=struct('name','WS1_','type',0,'data',[]);
Tag(length(Tag)+1)=struct('name','_DATASET_VERSION','type',5,'data',2);
Tag(length(Tag)+1)=struct('name','_DATASET_UID','type',10,...
  'data','2.16.124.113000.000000.00000000000000.000000000.00000000.0042');
Tag(length(Tag)+1)=struct('name','_DATASET_PATIENTUID','type',10,...
  'data','2.16.124.113000.000000.00000000000000.000000000.00000000.0011');
Tag(length(Tag)+1)=struct('name','_DATASET_PROCEDUREUID','type',10,...
  'data','2.16.124.113000.000000.00000000000000.000000000.00000000.0041');
Tag(length(Tag)+1)=struct('name','_DATASET_STATUS','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_DATASET_RPFILE','type',10,'data','default.rp');
Tag(length(Tag)+1)=struct('name','_DATASET_PROCSTEPTITLE','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_DATASET_PROCSTEPPROTOCOL','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_DATASET_PROCSTEPDESCRIPTION','type',10,...
  'data','');
Tag(length(Tag)+1)=struct('name','_DATASET_COLLECTIONDATETIME','type',10,...
  'data','Unknown');
Tag(length(Tag)+1)=struct('name','_DATASET_COLLECTIONSOFTWARE','type',10,...
  'data','Acq  ');
Tag(length(Tag)+1)=struct('name','_DATASET_CREATORDATETIME','type',10,...
  'data',sprintf('%d',floor(clock)));
Tag(length(Tag)+1)=struct('name','_DATASET_CREATORSOFTWARE','type',10,...
  'data','Acq ');
Tag(length(Tag)+1)=struct('name','_DATASET_KEYWORDS','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_DATASET_COMMENTS','type',10,...
  'data','Dummy infods.');
Tag(length(Tag)+1)=struct('name','_DATASET_OPERATORNAME','type',10,'data','');
Tag(length(Tag)+1)=struct('name','_DATASET_LASTMODIFIEDDATETIME','type',10,...
  'data',sprintf('%d',floor(clock)));
if exist_hc
  nominalPositions=0;   % Measured
else
  nominalPositions=1;   % Nominal
end
Tag(length(Tag)+1)=struct('name','_DATASET_NOMINALHCPOSITIONS','type',5,...
  'data',nominalPositions);
Tag(length(Tag)+1)=struct('name','_DATASET_COEFSFILENAME','type',10,...
  'data','ds.res4.scrr');
Tag(length(Tag)+1)=struct('name','_DATASET_SENSORSFILENAME','type',10,...
  'data','ds.res4.senres');
%Tag(length(Tag)+1)=struct('name','_DATASET_COEFSFILENAME','type',10,...
%  'data','/opt/ctf-5.1/hardware/M015/M015_1609.coef');
%Tag(length(Tag)+1)=struct('name','_DATASET_SENSORSFILENAME','type',10,...
%  'data','/opt/ctf-5.1/hardware/M015/M015_1609.sens');
Tag(length(Tag)+1)=struct('name','_DATASET_SYSTEM','type',10,'data','DSQ-2010');
Tag(length(Tag)+1)=struct('name','_DATASET_SYSTEMTYPE','type',10,'data','Untitled');
Tag(length(Tag)+1)=struct('name','_DATASET_LOWERBANDWIDTH','type',4,'data',0);
Tag(length(Tag)+1)=struct('name','_DATASET_UPPERBANDWIDTH','type',4,'data',...
  round(0.25*sample_rate));
Tag(length(Tag)+1)=struct('name','_DATASET_ISINDB','type',5,'data',0);

if exist_HLC
  HZ_MODE=5;
elseif exist_hc
  HZ_MODE=1;
else
  HZ_MODE=DATASET_HZ_UNKNOWN;
end
Tag(length(Tag)+1)=struct('name','_DATASET_HZ_MODE','type',5,'data',HZ_MODE);
Tag(length(Tag)+1)=struct('name','_DATASET_MOTIONTOLERANCE','type',4,'data',0.005);
Tag(length(Tag)+1)=struct('name','_DATASET_MAXHEADMOTION','type',4,'data',0.005);
Tag(length(Tag)+1)=struct('name','_DATASET_MAXHEADMOTIONTRIAL','type',7,'data',0);
Tag(length(Tag)+1)=struct('name','_DATASET_MAXHEADMOTIONCOIL','type',10,'data','1');
Tag(length(Tag)+1)=struct('name','EndOfParameters','type',-1,'data',[]);
Tag(length(Tag)+1)=struct('name','EndOfParameters','type',-1,'data',[]);
return
%  ************* End of function make_dummy_infods*********************************************
% **************************************************************************************


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Function combineMarkers  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrkNew=combineMarkers(mrkOld,nTrialOld,mrk,nTrialNew,res4);
%  Combines existing marker information with markers for new trial(s) being added to the
%  dataset ny addCTFtral.

%  Start by examining the new markers defined in structure mrk.
%  Examines structure ds to see if a sensible MarkerFile can be created from ds.mrk.
%  Output: mrkNew : New marker structure.  If the new markers are invlaid, or empty, then
%                   mrkNew=struct([]) is returned.

MarkerFileOK=~isempty(mrk);
if MarkerFileOK
  if isempty([mrk.trial]) | isempty([mrk.time])
    MarkerFileOK=0;  % Do not create MarkerFile.mrk if all of the marker classes are empty.
  else
    %  Are the markers appropriate?
    minMarkerTrial=min(min([mrk.trial]));
    maxMarkerTrial=max(max([mrk.trial]));
    minMarkerTime=min(min([mrk.time]));
    maxMarkerTime=max(max([mrk.time]));
    MarkerFileOK=(maxMarkerTrial<=nTrialNew & minMarkerTrial>=1 & ...
      maxMarkerTime<=(res4.no_samples/res4.sample_rate) & ...
      minMarkerTime>=(-res4.preTrigPts/res4.sample_rate));
    if ~MarkerFileOK
      fprintf(['addCTFtrial: Structure mrk cannot possibly be a set of markers ',...
          'for %d trial(s) in array(data).\n'],nTrialNew);
      fprintf(['          minMarkerTrial=%d (must be >=1)  ',...
          'maxMarkerTrial=%d (must be <=%d)\n'],...
        minMarkerTrial,maxMarkerTrial,ds.res4.no_trials);
      fprintf(['          minMarkerTime=%7.4f (must be >=%7.4f)  ',...
          'maxMarkerTrial=%7.4f  (must be <=%7.4f )\n'],...
        minMarkerTime,-res4.preTrigPts/res4.sample_rate,...
        maxMarkerTime,res4.no_samples/res4.sample_rate);
      fprintf('          MarkerFile.mrk will not be updated.\n');
    end
  end
end

if MarkerFileOK==0
  mrkNew=struct([]); % return an empty marker structure
  %  Check mrkOld to see if each existing trial had at least one marker.
  if isequal(unique([mrkOld.trial]),[1:nTrialOld])
    fprintf(['addCTFtrial: All of the trials in the existing data have at least one\n',...
        '             marked point,but the new trials have no marker points?\n']);
  end
  return
end

%  There are valid new markers.  Add the new markers to the existing marker structure,
%  and extend the marker definition if new markers have been defined.
%  Check the Name
mrkNew=mrkOld;
maxClassId=max([mrkNew.ClassId]);
for k=1:length(mrk)
  addq=1;
  %  Check if the marker is already defined.
  for q=1:length(mrkOld)
    if strcmp(mrkOld(q).Name,mrk(k).Name)
      mrkNew(q).trial=[mrkNew(q).trial mrk(k).trial+nTrialOld];
      mrkNew(q).time=[mrkNew(q).time mrk(k).time];
      addq=0;
      break;
    end
  end
  if addq   % Create a new marker class
    newClassId=~isfield(mrk(k),'ClassId');
    if ~newClassId
      newClassId=any(mrk(k).ClassId==[mrkNew.ClassId]);
    end
    if newClassId
      mrk(k).ClassId=max([mrkNew.ClassId])+1;
    end
    if ~isfield(mrk(k),'Color');mrk(k).Color='Red';end
    if ~isfield(mrk(k),'Comment');mrk(k).Comment=char([]);end
    if ~isfield(mrk(k),'ClassGroupId');mrk(k).ClassGroupId=3;end
    if mrk(k).ClassGroupId==0
      if ~isfield(mrk(k),'BitNumber') | ~isfield(mrk(k),'Polarity') | ...
          ~isfield(mrk(k),'Source') | ~isfield(mrk(k),'Threshold')
        fprintf(['addCTFtrial: Structure mrk defines ClassGroupId=0 marker %s, ',...
            'but not fields BitNumber, Polarity, Source or Threshold.\n'],mrk(k).Name);
        fprintf('             Marker set %s will be defined with ClassGroupId=3.\n',...
          mrk(k).Name);
        mrk(k).ClassGroupId=3;
      end
    end
    % ClassGroupId : make 4 fields empty.
    if mrk(k).ClassGroupId==3
      mrk(k).BitNumber=[];
      mrk(k).Polarity=[];
      mrk(k).Source=[];
      mrk(k).Threshold=[];
    end
    q0=length(mrkNew)+1;
    mrkNew(q0)=mrk(k);
    mrkNew(q0).trial=mrk(k).trial+nTrialOld;
    mrkNew(q0).time=mrk(k).time;
  end
  clear addq;
end
return
%%%%%%%%%%%%%% end of combineMarkers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function writeBadSegments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeBadSegments(badSegmentsFile,badSegments,nTrial,tTrial);

%  Creates a bad.segments file in a CTF data set from the information in structure
%  badSegments which is created by read_badSegments, or by the user.
%  If structure badSegments is empty, then the file is not created.

%  badSegments structure:
%    badSegments.trial = List of trial numbers   
%    badSegments.StartTime = List of bad segment start times (relative to trial).
%    badSegments.EndTime   = List of bad segment end times.

%  Check badSegmentsFile
if exist('badSegmentsFile')~=1;badSegmentsFile=char([]);end
if isempty(badSegmentsFile) | ~ischar(badSegmentsFile)
  fprintf('addCTFtrial(writeBadSegments):  Bad file name.\n');
  badSegmentsFile
  return
end

% Check that structure badSegments is defined correctly
if exist('badSegments')~=1 | exist('nTrial')~=1
  return
elseif ~isstruct(badSegments) | isempty(badSegments)
  return
elseif ~isfield(badSegments,'trial') | ~isfield(badSegments,'StartTime') | ...
    ~isfield(badSegments,'EndTime')
  return
elseif isempty(badSegments.trial) | isempty(badSegments.StartTime) | ...
    isempty(badSegments.EndTime)
  return
elseif ~isequal(size(badSegments.trial),size(badSegments.StartTime),...
    size(badSegments.EndTime))
  fprintf(['\naddCTFtrial (writeBadSegments): ',...
      'The fields of structure badSegments do not all have the same size.\n']);
  return
elseif any(badSegments.trial>nTrial) | any(badSegments.trial<1) | ...
    any(badSegments.StartTime<0) | any(badSegments.EndTime>tTrial)
  fprintf(['\naddCTFtrial (writeBadSegments): badSegments cannot possibly describe ',...
      'bad segments for these data.\n',...
      '\tmin(badSegments.trial)=%d max(badSegments.trial)=%d ',...
      'min(badSegments.StartTime)=%0.4f  max(badSegments.EndTime)=%0.4f s\n\t\tDataset:',...
      ' nTrial=%d  tTrial=%0.4f s\n'],...
    min(badSegments.trial),max(badSegments.trial),min(badSegments.StartTime),...
    max(badSegments.EndTime),nTrial,tTrial);
  fprintf('\t\tNew bad.segments file will not be created.\n\n');
  return
end

% Convert all fields to simple vectors
nSeg=prod(size(badSegments.trial));
trial=reshape(badSegments.trial,1,nSeg);
StartTime=reshape(badSegments.StartTime,1,nSeg);
EndTime=reshape(badSegments.EndTime,1,nSeg);

fid=fopen(badSegmentsFile,'w','ieee-be');
if fid<0
  fprintf('writeCTFds (writeBadSegments): Could not open file %s\n',badSegmentsFile);
  return
end
%  Extra tabs are inserted to reproduce the format of bad.segments files produced by
%  DataEditor (5.3.0-experimental-linux-20060918). 
fprintf(fid,'%0.6g\t\t%0.6g\t\t%0.6g\t\t\n',[trial;StartTime;EndTime]);
fclose(fid);
return
%%%%%%%%%  End of writeBadSegments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Function checkCls  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClassFileOK=checkCls(ds);

%  Examines structure ds to see if a sensible ClassFile can be created from ds.TrialClass.

ClassFileOK=isfield(ds,'TrialClass');
if ClassFileOK
  ClassFileOK=~isempty(ds.TrialClass);
end
if ClassFileOK
  %  Are the class trials appropriate?
  minClassTrial=[];
  maxClassTrial=[];
  for k=1:length(ds.TrialClass)
    maxClassTrial=max([maxClassTrial max(ds.TrialClass(k).trial)]);
    minClassTrial=min([minClassTrial min(ds.TrialClass(k).trial)]);
  end
  %  Create ClassFile.cls even when the trail classes are empty.
  if ~isempty(maxClassTrial)
    ClassFileOK=(maxClassTrial<=ds.res4.no_trials & minClassTrial>=1);
    if ~ClassFileOK
      fprintf(['\nwriteCTFds (checkCls): ds.TrialClass cannot possibly be a set of ',...
          'trial classes for array(data).\n          minClassTrial=%d (must be >=1)  ',...
          'maxClassTrial=%d (must be <=%d)\n'],...
        minClassTrial,maxClassTrial,ds.res4.no_trials);
      fprintf('          ClassFile.cls will not be created.\n');
    end
  end
end
return
%%%%%%%%%%%%%% end of checkCls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Function combineClasses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clsNew=combineClasses(clsOld,nTrialOld,cls,nTrialNew);
%  Combines existing trial class information with classes for new trial(s) that are
%  being added to the dataset by addCTFtral.

%  Start by examining the new markers defined in structure mrk.
%  Examines structure ds to see if a sensible CLassFile can be created from clsOld and cls.
%  Output: clsNew : New class structure.  If the new classes are invalid, or empty, then
%                   clsNew=struct([]) is returned.

ClassFileOK=~isempty(cls);
if ClassFileOK
  if isempty([cls.trial])
    ClassFileOK=0;  % Do not create ClassFile.cls if all of the trial classes are empty.
  else
    %  Are the markers appropriate?
    minClassTrial=min(min([cls.trial]));
    maxClassTrial=max(max([cls.trial]));
    ClassFileOK=(maxClassTrial<=nTrialNew & minClassTrial>=1);
    if ~ClassFileOK
      fprintf(['addCTFtrial: Structure cls cannot possibly be a set of classes ',...
          'for %d trial(s) in array data .\n'],nTrialNew);
      fprintf(['          minClassTrial=%d (must be >=1)  ',...
          'maxClassTrial=%d (must be <=%d)\n'],...
        minClassTrial,maxClassTrial,nTrialNew);
      fprintf('          ClassFile.cls will not be updated.\n');
    end
  end
end

if ClassFileOK==0
  clsNew=struct([]); % return an empty class structure
  %  Check clsOld to see if each existing trial had at least one class.
  if isequal(unique([clsOld.trial]),[1:nTrialOld])
    fprintf(['addCTFtrial: All of the trials in the existing data are assigned to at ',...
        'least one class\n',...
        '             but the new trials have no classes?\n']);
  end
  return
end

%  There are valid new classes.  Add the new classes to the existing class structure,
%  and extend the class definition if new classes have been defined.
clsNew=clsOld;
maxClassId=max([clsNew.ClassId]);
for k=1:length(cls)
  addq=1;
  %  Check if the class is already defined.
  for q=1:length(clsOld)
    if strcmp(clsOld(q).Name,cls(k).Name)
      clsNew(q).trial=[clsNew(q).trial cls(k).trial+nTrialOld];
      addq=0;
      break;
    end
  end
  if addq   % Create a new trial class
    newClassId=~isfield(cls(k),'ClassId');
    if ~newClassId  % Check if cls(k).ClassId is already in use.
      newClassId=any(cls(k).ClassId==[clsNew.ClassId]);
    end
    if newClassId
      cls(k).ClassId=max([clsNew.ClassId])+1;
    end
    if ~isfield(cls(k),'Color');cls(k).Color='Red';end
    if ~isfield(cls(k),'Comment');cls(k).Comment=char([]);end
    if ~isfield(cls(k),'ClassGroupId');cls(k).ClassGroupId=3;end
    q0=length(clsNew)+1;
    clsNew(q0)=cls(k);
    clsNew(q0).trial=cls(k).trial+nTrialOld;
  end
  clear addq;
end
return
%%%%%%%%%%%%%% end of combineClasses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Function writeClassFile    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeClassFile(ClassFile,TrialClass);

%   Write the ClassFile of a CTF data set.
%   The CLassFile allows a user to store a list of trial classifications in a data set.
%   The ClassFile format is defined in document CTF MEG File Formats, PN900-0088.
%   This format is rigid. 

%   Inputs :
%      ClassFile : marker file including the full path and extension .mrk.
%      TrialClass : Structure creted by read_ClassFile.

%  Output : ClassFile.cls.

%  Check input TrialClass.
if exist('TrialClass')~=1;TrialClass=[];end
if isempty(TrialClass) | ~isstruct(TrialClass)
  fprintf('writeCTFds (writeClassFile):  TrialClass is empty or is not a structure.\n');
  TrialClass
  return
end

%  Check ClassFile
if exist('ClassFile')~=1;ClassFile=char([]);end
if isempty(ClassFile) | ~ischar(ClassFile)
  fprintf('writeCTFds (writeClassFile):  Bad file name.\n');
  ClassFile
end

fid=fopen(ClassFile,'w','ieee-be');
if fid<0
  fprintf('writeCTFds (writeClassFile): Could not open file %s\n',ClassFile);
  return
end

nClass=length(TrialClass);

%  Generate datasetname from ClassFIle.
ksep=max([0 strfind(ClassFile,filesep)]);
datasetname=ClassFile(1:ksep-1);
if isempty(datasetname);datasetname=cd;end

fprintf(fid,'PATH OF DATASET:\n%s\n\n\n',datasetname);
fprintf(fid,'NUMBER OF CLASSES:\n%d\n\n\n',nClass);

for k=1:nClass
  if k==1           % Add sign character to make output match the output of Acq.
    sgn=char([]);   % There should be no real significance to this.
  else              % Why does DataEditor places the + sign only on ClassID 2,3,...?
    sgn='+';
  end
  No_of_Trials=prod(size(TrialClass(k).trial));
  fprintf(fid,'CLASSGROUPID:\n%s%d\nNAME:\n%s\nCOMMENT:\n%s\n',...
    sgn,TrialClass(k).ClassGroupId,TrialClass(k).Name,TrialClass(k).Comment);
  fprintf(fid,'COLOR:\n%s\nEDITABLE:\n%s\nCLASSID:\n%s%d\nNUMBER OF TRIALS:\n%d\n',...
    TrialClass(k).Color,TrialClass(k).Editable,sgn,TrialClass(k).ClassId,No_of_Trials);
  fprintf(fid,'LIST OF TRIALS:\nTRIAL NUMBER\n');
  % Subtract one so trial numbering starts at 0 in ClassFile.cls
  fprintf(fid,'%20d\n',reshape(TrialClass(k).trial,1,No_of_Trials)-1); 
  fprintf(fid,'\n\n');
end
fclose(fid);
return
%%%%%%%%%%%%%% end of writeClassFile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Function updateCreatorSoftware   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ds=updateCreatorSoftware(ds,creatorSoftware,originalCreatorSoftware);

%  Changes the creator software fields of the ds structure, and returns ds=struct([]) if 
%  field infods is missing from structure ds, or if it thinks that the dataset was not
%  created by the MATLAB code specified by originalCreatorSoftware (i.e. writeCTFds).

%  Makes sure that creatorSoftware info is changed in structure infods, newds and res4.

%  Inputs: ds : ds structure produced by readCTFds.
%          creatorSoftware : Character string indicating that the data set was
%              modified by addCTFtrial.  Added to infods tags listed in addCreatorTag and
%              appName field of res4.
%          originalCreatorSoftware : Character string indicating that the data set was
%              created by writeCTFds.  Added to infods tags listed in addCreatorTag and
%              appName field of res4.

%  Output: ds : ds structure with changes to ds.infids, ds.res4, ds.newds.  
%               If the input ds structure does not indicate that iot was created by
%               writeCTFds (originalCreatorSoftware), ds=struct([]) is returned.

%  Adds comment (clinical use message) and creator software name to infods, newds and hist files.

if ~isfield(ds,'infods');
  fprintf(['addCTFtrial (updateCreatorSoftware): Structure ds does not include',...
      ' field infods.\n']);
  whos ds
  ds=struct([]);  % Force an error in the calling program
  return
end

%  The infods tags and red4 fields that display creator software must also indicate that
%  the dataset was created by originalCreatorSoftware.

%  infods tags that will display the creator software.
addCreatorTag=strvcat('_PROCEDURE_COMMENTS','_DATASET_STATUS',...
  '_DATASET_COLLECTIONSOFTWARE','_DATASET_CREATORSOFTWARE');
%  res4 text fields that will display the creator software.
addCreatorField=strvcat('appName','run_description');
addCreatorLength=[256 -1]'; % -1 indicates variable length

if exist('creatorSoftware')~=1;
  creatorSoftware=char([]);
elseif ~ischar(creatorSoftware)
  creatorSoftware=char([]);
else
  creatorSoftware=deblank(creatorSoftware);
end

if exist('originalCreatorSoftware')~=1;
  originalCreatorSoftware=char([]);
elseif ~ischar(originalCreatorSoftware)
  originalCreatorSoftware=char([]);
else
  originalCreatorSoftware=deblank(originalCreatorSoftware);
end
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that the string originalCreatorSoftware appears in the infods and res4 files.

tagName=getArrayField(ds.infods,'name')';
if length(ds.infods)==1;tagName=tagName';end

%  Update the infods fields.
for k=1:size(addCreatorTag,1)
  q=strmatch(deblank(addCreatorTag(k,:)),tagName);
  if isempty(q);
    fprintf('addCTFtrial (updateCreatorSoftware): Tag %s is missing from ds.infods.\n',...
      deblank(addCreatorTag(k,:)));
    ds=struct([]);
    return
  elseif length(q)>1
    fprintf(['addCTFtrial (updateCreatorSoftware): Tag %s appears %d times in ',...
        'ds.infods.\n'],deblank(addCreatorTag(k,:)),length(q));
    ds=struct([]);
    return
  else
    strng=ds.infods(q).data;
    % Strng must contain originalCreatorSoftware
    nChar=length(strng);
    if length(originalCreatorSoftware)>0 
      qpos=strfind(strng,originalCreatorSoftware);
      if isempty(qpos)
        fprintf(['addCTFtrial (updateCreatorSoftware): String %s does not appear in',...
            ' ds.infods tag %s\n'],originalCreatorSoftware,addCreatorTag(k,:));
        fprintf(['                                     This dataset was not created ',...
            'by CTF MATLAB Export software.\n']);
        ds=struct([]);
        return
      end
      qpos=qpos+length(originalCreatorSoftware)-1;
    else
      qpos=length(strng);
    end
    %  Add creatorSoftware string if it is not already present.
    if isempty(strfind(strng,creatorSoftware))
      ds.infods(q).data=[strng(1:qpos) ' ' creatorSoftware ' ' strng(qpos+1:nChar)];
    end
  end
end

%  Update the res4 fields: add creatorSoftware message.  
%  Don't check field length.  Truncate later.
for q=1:size(addCreatorField,1);
  strng=deblank(getfield(ds.res4,deblank(addCreatorField(q,:))));
  nChar=length(strng);
  if length(originalCreatorSoftware)>0 
    qpos=strfind(strng,originalCreatorSoftware);
    if isempty(qpos)
      fprintf(['addCTFtrial (updateCreatorSoftware): String %s does not appear in',...
          ' ds.res4.%s\n'],originalCreatorSoftware,deblank(addCreatorField(q,:)));
      fprintf(['                                      This dataset was not created ',...
          'by CTF MATLAB Export software.\n']);
      ds=struct([]);
      return
    end
    qpos=qpos+length(originalCreatorSoftware);
  else
    qpos=nChar;
  end
  %  Add creatorSoftware string if it is not already present.
  if isempty(strfind(strng,creatorSoftware))
    newStrng=[strng(1:qpos) ' ' creatorSoftware ' ' strng(qpos+1:nChar)];
    ds.res4=setfield(ds.res4,deblank(addCreatorField(q,:)),newStrng);
  end
end
clear q strng nChar strng0 ns newStrng;

% Update res4.run_description
ds.res4.run_description=[deblank(ds.res4.run_description),char(0)];
ds.res4.rdlen=length(ds.res4.run_description);

%  Truncate the .res4. fields.  Leave room for a final char(0).
for q=1:size(addCreatorField,1);
  strng=deblank(getfield(ds.res4,deblank(addCreatorField(q,:))));
  if length(strng)>addCreatorLength(q)-1 & addCreatorLength(q)>0
    ds.res4=setfield(ds.res4,deblank(addCreatorField(q,:)),...
      strng(length(strng)+[-addCreatorLength(q)+2:0]));
  end
end
clear q strng;

%  Add the creator and comment information to the newds file
if isfield(ds,'newds');
  nChar=length(ds.newds);
  %  Keyword positions
  aNpos=max(strfind(ds.newds,[char(9) 'appName:' char(9)]));
  if isempty(aNpos)
    fprintf(['addCTFtrial (updateCreatorSoftware): Keyword ''appName'' ',...
        'does not appear in ds.newds.\n',...
        '         set ds.newds=[].\n']);
    ds.newds=char([]);
  else
    eol=max([11 min(strfind(ds.newds(aNpos+1:length(ds.newds)),char(10)))]);
    wPos=min(strfind(ds.newds(aNpos+[1:eol]),originalCreatorSoftware));
    addPos=min(strfind(ds.newds(aNpos+[1:eol]),creatorSoftware));
    if isempty(addPos)
      if isempty(wPos)
        ds.newds=[ds.newds(1:(aNpos+eol-1)) ' ' creatorSoftware ...
            ds.newds((aNpos+eol):length(ds.newds))];
      else
        ds.newds=[ds.newds(1:(aNpos+wPos+length(originalCreatorSoftware)-1)) ', ' ...
            creatorSoftware ...
          ds.newds((aNpos+wPos+length(originalCreatorSoftware)):length(ds.newds))];
      end
    end
    clear eol wPos addPos;
  end
  clear nChar aNpos;
end

return
%%%%%%%%%%%%%% End of update_creatorSoftware %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Function updateHLC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [infods,status]=updateHLC(infods,HLCdata,HLC0,nTrial0);
%  Adds notes and  continuous-head-localization tags to .infods file.
%  Updates MAXHEADMOTION tags when called from addCTFtrial

%  Inputs : infods : Structure of tags to be written to output infods file.
%           HLCdata : 0 or [] or not defined - There are no HLC channels
%                     (head-coil position) in the data.  Just pass the existing tags to
%                     the output dataset.
%           HLCdata : real array : size(HLCdata)=[samplesPerTrial 3 Ncoil nTrial]
%                     Calculate the head motion tags.
%           HLC0 : Head coil positions (m) in dewar coordinates at the start of the run.
%                  size(HLC0)=[3 Ncoil]
%           nTrial0=No. of trials of data already written to the dataset.


%  Outputs : status :  0: Everything looks OK
%                      1: Couldn't find _DATASET tags.
%                     <0: There is a problem with the new infods file.

%  Calls: 
%     - getArrayField: Extracts field names from structure array S.

% Defaults for the head localization tags.
% Might do this better with a structure.

HZ_MODE_UNKNOWN=2^31-1;   % Bug in DataEditor and acquisition software, Oct. 2006
%                           Should be -1.

%  New dataset tags for HLC specification.
addTag=strvcat('_DATASET_HZ_MODE',...
  '_DATASET_MOTIONTOLERANCE',...
  '_DATASET_MAXHEADMOTION',...
  '_DATASET_MAXHEADMOTIONTRIAL',...
  '_DATASET_MAXHEADMOTIONCOIL');
addTagType=[5 4 4 7 10]';

if exist('nTrial0')~=1
  nTrial0=0;
elseif isempty(nTrial0)
  nTrial0=0;
elseif min(nTrial0)<0
  nTrial0=0;
end
%  Check HLCdata array.  HLCdata=[] indicates no continuous head localization.
if exist('HLCdata')~=1
  HLCdata=[];
elseif ~isempty(HLCdata)
  [samplesPerTrial nDim Ncoil nTrial]=size(HLCdata);
  if nDim~=3;
    fprintf('addCTFtrial (updateHLC): size(HLCdata)=[');fprintf(' %d',size(HLCdata));
    fprintf(']\n');
    fprintf('             nDim=%d??   Setting HLCdata=[].\n',nDim);
    HLCdata=[];
  end
end
if exist('HLC0')~=1
  HLC0=[];
elseif ~isequal(size(HLC0),[3 size(HLCdata,3)])
  HLC0=[];
end
if isempty(HLC0) & ~isempty(HLCdata)
  HLC0=squeeze(HLCdata(1,:,:,1));
end

%  Assume that all the tag names are different.  There is no checking when a tag listed in
%  addTag matches more than one entry in array name.
name=getArrayField(infods,'name')';  %size(name)=[nTag lengthTag]
DATASET_tags=strmatch('_DATASET',name)';  %  size(DATASET_tags)=[1 nTag]
if isempty(DATASET_tags)
  status=1;  % No _DATASET tags.  Don't add anything.
else
  status=0;
  if isempty(HLCdata)
    TagValue(1)=HZ_MODE_UNKNOWN;
    TextTagValue=char(0);
    addTag=addTag(1,:);   % Remove the other HLC tags
    addTagType=addTagType(1);
  else % ~isempty(HLCdata)
    %  Remove HLC offsets.
    for q=1:Ncoil
      for k=1:3
        HLCdata(:,k,q,:)=HLCdata(:,k,q,:)-HLC0(k,q);  
      end
    end
    %  Calculate motions as displacement from the start of the dataset.
    absMotion=squeeze(sqrt(sum(HLCdata.^2,2)));  %size(absMotion)=[samplesPerTrial Ncoil nTrial]
    maxCoilMotion=squeeze(max(absMotion,[],1));  % size(maxCoilMovement)=[Ncoil nTrial]
    maxHeadMotion=max(max(maxCoilMotion));
    [mx maxHeadMotionCoil]=max(max(maxCoilMotion,[],2));
    [mx maxHeadMotionTrial]=max(max(maxCoilMotion,[],1));
    %  Create a list of head motion tag values
    TagValue(1)=5;      % Indicates continuous head localization
    TagValue(2)=max(2*maxHeadMotion,0.02); % _DATASET_MOTIONTOLERANCE
    TagValue(3)=maxHeadMotion;
    TagValue(4)=maxHeadMotionTrial-1+nTrial0;  % _DATASET_MAXHEADMOTIONTRIAL
    TextTagValue=strvcat(char(zeros(4,1)),sprintf('%d',maxHeadMotionCoil));
    
    %  Does the existing infods have the head-motion tags?
    TagNo2=strmatch(deblank(addTag(2,:)),name,'exact');
    TagNo3=strmatch(deblank(addTag(3,:)),name,'exact');
    TagNo4=strmatch(deblank(addTag(4,:)),name,'exact');
    TagNo5=strmatch(deblank(addTag(5,:)),name,'exact');
    if ~isempty(TagNo2) & ~isempty(TagNo3) & ~isempty(TagNo4) & ...
        ~isempty(TagNo5) & nTrial0>0
      %  Structure infods already has the head-motions tags from previous data.
      TagValue(2)=max(TagValue(2),infods(TagNo2).data);
      if TagValue(3)<=infods(TagNo3).data
        TagValue(3)=infods(TagNo3).data;
        TagValue(4)=infods(TagNo4).data;
        TextTagValue=strvcat(char(zeros(4,1)),infods(TagNo5).data);
      end
    elseif (~isempty(TagNo2) | ~isempty(TagNo3) | ~isempty(TagNo4) | ...
        ~isempty(TagNo5)) & nTrial0>0
      fprintf(['addCTFtrial (updateHLC):  Some, but not all, of the ds.infods CHL ',...
          'tags are defined.\n',...
          '            Infods will be defines with information from the %d new ',...
          'trial(s) being added.\n'],size(HLCdata,4));
    end
    TagValue=[TagValue 0]; % Placeholder only since the 5th tag is a text string.
  end
  
  % Add or insert tags.
  for q=1:size(addTag,1)
    nTag=length(infods);
    tagName=deblank(addTag(q,:));
    TagNo=strmatch(tagName,name,'exact')';
    if isempty(TagNo)  % Insert a tag at the end of the _DATASET tags.
      TagNo=max(DATASET_tags)+1;
      infods((TagNo+1):(nTag+1))=infods(TagNo:nTag);
      name=strvcat(name(1:TagNo-1,:),tagName,name(TagNo:nTag,:));
      DATASET_tags=[DATASET_tags TagNo];
    end
    if addTagType(q)~=10
      infods(TagNo)=struct('name',tagName,'type',addTagType(q),'data',TagValue(q));
    else
      infods(TagNo)=struct('name',tagName,'type',addTagType(q),...
        'data',deblank(TextTagValue(q,:)));
    end
  end  % End loop over head position and head motion tags.
  clear q TagNo TextTagValue TagValue;
end  % End section add _DATASET tags

return
%%%%%%%%%%%%%% End of updateHLC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Function getArrayField  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=getArrayField(S,yfield)

%  Extracts one field of a structure array.
%  Inputs:  S : structure array.
%           yfield: name of field of S (string)

%  Output: x.  x(:,n)=S(n).yfield;  size(x)=[size(yfield) length(S)]

%  Field sizes: Two options:
%    1. size(S(n).yfield) is the same for all n.  Any size allowed.
%    2. S(n).yfield is a 2-D array for all structure elements S(n).
%      Then yfield of different array elements S(n) can have different sizes.
%      The array returned will be large enough to accomodate all of the data.

sizeS=size(S);
nS=prod(sizeS);
S=reshape(S,1,nS);

%  Determine which array-size option to use.
sizey=size(getfield(S,{1},yfield));
option1=1;
option2=(length(sizey)==2);
for n=2:nS
  sizeyn=size(getfield(S,{n},yfield));
  option1=option1 & isequal(sizey,sizeyn);
  option2=option2 & length(sizeyn)==2;
  if option2
    sizey=max([sizey;sizeyn],[],1);
  end
end

if option1  % All fields have the same size
  nY=prod(sizey);
  if isnumeric(getfield(S,{1},yfield))
    x=zeros(nY,nS);
  elseif ischar(getfield(S,{1},yfield))
    x=char(zeros(nY,nS));
  end
  for n=1:nS
    x(:,n)=reshape(getfield(S,{n},yfield),nY,1);
  end
  x=reshape(x,[sizey nS]);
elseif option2  % All fields have only two dimensions
  if isnumeric(getfield(S,{1},yfield))
    x=zeros([sizey nS]);
  elseif ischar(getfield(S,{1},yfield))
    x=char(zeros([sizey nS]));
  end
  for n=1:nS
    y=getfield(S,{n},yfield);
    sizeyn=size(y);
    x(1:sizeyn(1),1:sizeyn(2),n)=y;
  end
else
  fprintf(['getArrayField: Field %s of the structure array has >2 dimensions ',...
      'and not all field have the same size.\n'],yfield);
  x=[];
end
x=squeeze(x);
return
%%%%%%%%  End of getArrayField  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
