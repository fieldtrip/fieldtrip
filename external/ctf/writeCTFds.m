function ds=writeCTFds(datasetname,ds,data,unit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                %
%   This program creates datasets that can be analyzed by CTF software.          %
%                                                                                %
%   Datasets created by this program MUST NOT BE USED FOR CLINICAL APPLICATIONS. %
%                                                                                %
%   Please do not redistribute it without permission from VSM MedTech Ltd.       %
%                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Author : Harold Wilson

%  Version 1.3  5 October 2007  Spelling errors in some variables corrected.
%  Version 1.2  24 April 2007  Modified to write both MEG and fMEG .hc files.

%  writeCTFds.m    Version 1.1  Prepares a CTF-format data set. MATLAB to CTF conversion.
%                  Adds record of filter changes to hist and newds files.
%                  Adds no-clinical-use messages.
%                  Creates multiple meg4 files when the total size of the data array
%                  >536870910 elements.

%  Operation :
%    1.  User reads a data set using readCTFds and getTrial2.
%        size(data)= [SAMPLES, CHANNELS, TRIALS].
%    2.  After working on data in MATLAB, user adjusts ds structure to reflect changes.
%        (writeCTFds will adjust the number of channels, channel names and trial structure.)
%    3.  This program then creates a new CTF data set.
%        datasetname must include the complete path description.
%        If a data set with name datasetname already exists, writeCTFds will issue an error message.
%        The new directory contains files for each field of structure ds.  If the field is
%        missing no file is created.  If the field is empty, an empty file is created.
%        Files default.* are not created by writeCTFds.

%    4.  The following fields of structure ds.res4 are modified based on the size of array data :
%             no_samples 
%             no_channels
%             no_trials.
%        If size(ds.res4.chanNames,1)<no_channels, additional channel names are created
%        as required.  The additional channels are called MXT%%%.  It is assumed that 
%        there will be <1000 new channels.

%  Inputs : datasetname : Output dataset including path.  Extension '.ds' is optional.
%           ds : Structure produced by readCTFds.
%           data : MEG data array.  size(data)=[no_samples no_channels no_trials]
%                  Array data may be single or double.
%           unit :  Determines the unit of the SQUID and EEG signals:
%                   If unit is missing or unit==[], the unit is set to 'fT'.
%                   'ft' or 'fT' : Convert to fT (MEG), uV (EEG)
%                   't'  or 'T'  : Convert to  T (MEG),  V (EEG)
%                   'phi0' : Convert to phi0 (MEG), uV (EEG)
%                   'int': Read plain integers from *.meg4-file

%  Outputs : - ds  :  The ds structure of the output data set.
%            - datasetout : the name of the output data set.
%            - A data set.  The .hz and .hz2 subdirectories are not included.

%  Function calls
%   Included in this listing:
%      - check_senres: Does simple checks on the fields of the senres table.
%      - writeHc:    Creates the new .hc file
%      - checkMrk:   Checks structure ds.mrk to see if it is valid marker set.
%      - writeEEG:     Creates the new .EEG file.
%      - writeBadSegments:  Creates the new bad.segments file.
%      - writeClassFile:    Creates the new ClassFile.cls file.
%      - writeVirtualChannels:  Creates the new VirtualChannels file.
%      - updateDescriptors:  Adds non-clinical-use and creation software messages to
%                 infods, res4, newds and hist fields of ds.
%      - updateHLC: Adds head-coil movement information to infods.
%      - updateDateTime : Resets dattime fields of res4 and infods.
%      - updateBandwidth: Resets bandwidth of newds and infods.  Adds res4 filter 
%                          description to ds.hist.
%      - getArrayField : Extracts one field of a structure array to make it easier to
%                        manipulate.
%      - writeMarkerFile:  Creates the new MarkerFile.mrk file.
%      - writeCPersist:     Creates the new .acq and .infods files.

%   Other calls:
%      - writeRes4: Writes the new .res4 file.

%  Output files are opened with 'w' permission and 'ieee-be' machine format in order
%  to be compatible with the Linux acquisition and analysis software.  Do not open files
%  with 'wt' permission because this will add an extra char(13) byte at the end of each
%  line of text.

persistent printWarning bandwidthMessage
delim=filesep;

if nargin==0 & nargout==0  %  Print a version number
  fprintf(['\twriteCTFds: Version 1.3   5 October 2007   ',...
      'Creates v4.1 and v4.2 CTF data sets.\n',...
      '\tCall: ds=writeCTFds(datasetname,ds,data,unit);\n',...
      '\t\tdatasetname = Name of the new dataset including the path.\n',...
      '\t\tds = Structure describing the new dataset (ds.hc must be v1.2 format).\n',...
      '\t\tdata = data that will be written to the new dataset .meg4 file.\n',...
      '\t\tunit = physical units of the data.\n\n']);
  return
end

%  Allowed 8-byte headers for res4 and meg4 files.
res4_headers=strvcat(['MEG41RS',char(0)],['MEG42RS',char(0)]);
meg4_headers=strvcat(['MEG41CP',char(0)],['MEG42CP',char(0)]);

maxMEG4Size=2^31;  % Maximum MEG$ file in bytes.  (Limit set by Linux software)

MAX_COILS=8;      % Parameter that determines the size of the ds.res.senres structure.
lenSensorName=32;  % Channel names must be 32 characters
printDefaultBandwidthMessage=0;  % Print a message about default bandwidth?
default_flp=0.25;  % Default is flp=0.25*ds.res4.sample_rate
clinicalUseMessage='NOT FOR CLINICAL USE';
creatorSoftware='writeCTFds';  % Added to .infods file
meg4ChunkSize=2^20; % Write new .meg4 file in chunks of 4*meg4ChunkSize bytes.
DATASET_HZ_UNKNOWN=round(2^31-1);  % Peculiar requirement of DataEditor as of 23 Oct. 2006

if ~exist('clinicalUseMessage');
  clinicalUseMessage=char([]);
end

%  Check inputs
if nargin<3
  fprintf(['\nwriteCTFds: Must supply inputs datasetname,ds,data.  ',...
      'Only %d input arguments are present.\n\n'],nargin);
  ds=-1;  % Force an error in the calling program.
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
    fprintf(['\nwriteCTFds : unit=%s    Not a valid option.  Must be ',...
        '''fT'', ''T'', ''phi0'' or ''int''\n\n'],unit);
    ds=-1;  % Force an error in the calling program.
    return
  end
end

% Check argument type
if ~isstruct(ds) | ~isnumeric(data) | ~ischar(unit) | ~ischar(datasetname)
  fprintf('\nwriteCTFds: Some of the inputs are the wrong type.\n');
  whos datasetname ds data unit;
  ds=-1;
  return
elseif ~isfield(ds,'res4') | ~isfield(ds,'meg4')
  fprintf('\nwriteCTFds: Fields res4 and meg4 must be present in structure ds.\n\n');
  ds  % List the fields of structure ds.
  ds=-1;  % Force an error in the calling program.
  return
end

% Refuse to write a data set with balanced reference gradiometers.
balancedGref=0;
for k=1:ds.res4.no_channels
  balancedGref=(ds.res4.senres(k).sensorTypeIndex==1 & ds.res4.senres(k).grad_order_no~=0);
end
if balancedGref
  fprintf('\nwriteCTFds: ds.res4.senres indicates balanced reference gradiometers.\n\n');
  ds=-1;  % Force an error in the calling program.
  return
end
clear k balancedGref;

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

%  Save the name already in structure ds, and change to the new datset name.
if isfield(ds,'baseName')
  olddatasetname=[ds.baseName,'.ds'];
  if isfield(ds, 'path')
    olddatasetname=[ds.path,olddatasetname];
  end
else
  olddatasetname='  None  ';
end
ds.path=path;
ds.baseName=baseName;

% Does the dataset already exist?
if exist(datasetname)==7
  fprintf('\nwriteCTFds: Dataset %s already exists.  Use a different name.\n\n',...
    datasetname);
  ds=-1;  % Force an error in the calling program.
  return
end

if size(ds.res4.chanNames,2)~=lenSensorName
  fprintf(['\nwriteCTFds : size(ds.res4.chanNames)=[%d %d] ?  Must ',...
      'have %d-character channel names.\n\n'],size(ds.res4.chanNames),lenSensorName);
  ds=-1;
  return
end

%  Check that the channel names have a sensor-file identification extensions.
%  If it is missing, print a warning message.
%  Sensor type indices : SQUIDs (0:7), ADCs (10), DACs(14), Clock (17), HLC (13,28,29)
%  See Document CTF MEG File Formats (PN 900-0088), RES4 File Format/
for index=[0:7 10 13 14 17 28 29]
  for k=find([ds.res4.senres.sensorTypeIndex]==index);
    if isempty(strfind(ds.res4.chanNames(k,:),'-'))
      fprintf(['writeCTFds: Channel %3d  %s     No sensor-file identification.',...
          ' (''-xxxx'' appended to channel name).\n',...
          '\t\tSome CTF software may not work with these channel names.\n'],...
        k,deblank(ds.res4.chanNames(k,:)));
      break;
    end
  end
  if isempty(k)||isempty(strfind(ds.res4.chanNames(k,:),'-'));break;end
end
clear index k chanName;

%  Update the data description in the ds.res4 structure.
[nSample, nChan, trials]=size(data);

%  Update ds.res4 fields to match the size of array data.
ds.res4.no_trials=trials;
ds.res4.no_channels=nChan;
ds.res4.no_samples=nSample;
ds.res4.epoch_time=nSample*trials/ds.res4.sample_rate;

%  Check if channels have been added or removed from array data.
[no_chanNames len_chanNames]=size(ds.res4.chanNames);
if no_chanNames<nChan
  % Assume that new channel are at the end of the data set.  Add a fake extension.
  for kx=1:(nChan-no_chanNames)
    ds.res4.chanNames=...
      strvcat(ds.res4.chanNames,['MXT' num2str(kx,'%3.3d') '-0001' char(0)]);
  end
  fprintf('\tAdded %d SQUID channels to the end of ds.res4.chanNames table.\n',...
    nChan-no_chanNames);
elseif no_chanNames>nChan
  fprintf(['\nlength(chanNames)=%d, but only %d channels of data.  ',...
      'writeCTFds cannot tell which channel names to remove.\n\n'],no_chanNames,nChan);
  ds=-1;
  return
end
clear no_chanNames len_chanNames;

% The senres table may have been changed, especially if channels are added or removed.
% Check structure senres, print error messages, and possibly fix the errors.
[ds.res4.senres,status]=check_senres(ds.res4.senres,ds.res4.no_channels);
if status<0;ds=-1;return;end
clear status;

%  Check that ds.res4.numcoef is the size of structure array ds.res4.scrr.
if isfield(ds.res4,'scrr')
  if ~isequal(size(ds.res4.scrr,2),ds.res4.numcoef)
    fprintf('Error in ds.res4:  ds.res4.numcoef=%d, but size(ds.res4.scrr)=[',...
      ds.res4.numcoef);
    fprintf(' %d',size(ds.res4.scrr));fprintf(']  ?\n');
    return
  end
elseif ds.res4.numcoef~=0
  fprintf(['Error in ds.res4:  ds.res4.numcoef=%d,',...
      ' but scrr is not a field of ds.res4\n'],ds.res4.numcoef);
  return
end

%  Before converting data to integers, save HLC channels for motion analysis in function
%  make_new_infods.  Pass HLCdata to function make_new_infods
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
  HLCdata=reshape(double(data(:,HLClist,:)),ds.res4.no_samples,3,coil,ds.res4.no_trials);
  clear coil k HLClist;
end

%  Convert data to integers because CTF data sets are stored as raw numbers and
%  not as physical qunatities.  The res4 file contains the calibrations for 
%  converting back to physical units.  Array data may be single precision, so 
%  convert to double before doing any adjustments to the data.
%  Already checked that unit is valid.

if strcmp(unit,'int')
  data=reshape(data,nSample,nChan*trials);
  for k=1:nChan*trials
    if strcmp(class(data),'single')
      data(:,k)=single(round(double(data(:,k))));
    else
      data(:,k)=round(double(data(:,k)));
    end
  end
  data=reshape(data,nSample,nChan,trials);
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
    for kt=1:trials
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

%  Create the output dataset
[status,msg]=mkdir(path,[baseName '.ds']);
if status==0 | ~isempty(msg)
  fprintf('\nwriteCTFds: Failed to create directory.\n');
  fprintf('            [status,msg]=mkdir(%s)\n',datasetname);
  fprintf('            returns status=%d, msg=%s',status,msg);fprintf('\n\n');
  ds=-1;
  return
end
clear msg status;


%  Write the data file (meg4 file).  Check ds.meg4.header and make sure that 
%  the output .meg4 file has an acceptable header.
headerMessage=1;
if ~isfield(ds.meg4,'header');  %  If ds.meg4.header is missing, add it.
  nChar=length(deblank(meg4_headers(1,:)));
  ds.meg4.header=[meg4_headers(1,1:min(7,nChar)) char(zeros(1,7-nChar))];
elseif isempty(strmatch(ds.meg4.header(1:7),meg4_headers(:,1:7),'exact'))
  ds.meg4.header=meg4_headers(1,1:7);
else
  headerMessage=0;
end
if headerMessage;
  fprintf('writeCTFds: Set ds.meg4.header=%s\n',ds.meg4.header);
end
clear headerMessage;

if isempty(printWarning)
  fprintf(['\nwriteCTFds: The data you are writing have been processed by software not\n',...
      '\tmanufactured by VSM MedTech Ltd. and that has not received marketing clearance\n',...
      '\tfor clinical applications.  These data should not be later employed for clinical\n',...
      '\tand/or diagnostic purposes.\n\n']);
  printWarning=1;
end

%  Write the meg4 file(s).  If there are more than maxMEG4Size-8 bytes, then additional meg4
%  files will be created.
%  Convert data to a 1-D array
ndata=prod(size(data));
data=reshape(data,ndata,1);
ptsPerTrial=nSample*nChan;
maxPtsPerFile=ptsPerTrial*floor((maxMEG4Size-8)/(4*ptsPerTrial));
pt=0;   % Last point written to the output file(s).
while pt<ndata
  endPt=pt+min(ndata-pt,maxPtsPerFile);
  if pt==0
    meg4Ext='.meg4';
  else
    meg4Ext=['.',int2str(floor(pt/maxPtsPerFile)),'_meg4'];
  end
  fidMeg4=fopen([path,baseName,'.ds',delim,baseName,meg4Ext],'w','ieee-be');
  fwrite(fidMeg4,[ds.meg4.header(1:7),char(0)],'uint8');
  while pt<endPt
    pt1=min(pt+meg4ChunkSize,endPt);                   %  Convert to double in case data is
    fwrite(fidMeg4,double(data((pt+1):pt1)),'int32');  %  is single and write in short 
    pt=pt1;                                            %  pieces. 
  end
  fclose(fidMeg4);
end
%  Update the .meg4 part of structure ds.
ds.meg4.fileSize=4*ndata+8*(1+floor(ndata/maxPtsPerFile));
clear data pt pt1 ndata fidMeg4 ptsPerTrial maxPtsPerFile meg4Ext;

%  Add dataset names to .hist
if ~isfield(ds,'hist');ds.hist=char([]);end
ds.hist=[ds.hist char(10) char(10) datestr(now) ' :' char(10) ...
    '      Read into MATLAB as data set ' olddatasetname char(10) ...
    '      Rewritten by writeCTFds as data set ' datasetname char(10)];

%  If infods doesn't exist or is empty create it.
if ~isfield(ds,'infods');ds.infods=[];end
if isempty(ds.infods)
  ds.infods=make_dummy_infods(isfield(ds,'hc'),~isempty(HLCdata),ds.res4.sample_rate);
end

%  Update text fields of res4,infods, newds and hist.
ds=updateDescriptors(ds,clinicalUseMessage,creatorSoftware);
%  Add HLC data to infods
ds.infods=updateHLC(ds.infods,HLCdata);

%  Analyze structure array ds.res4.filters to make text info for .hist file and
%  bandwidth parameters for .newds file.
fhp=0;  % High-pass cutoff assuming no filters
flp=default_flp*ds.res4.sample_rate;   % Assumed lowpass cutoff.  SHOULD THIS BE CHANGED?
if isempty(bandwidthMessage) & printDefaultBandwidthMessage
  fprintf('writeCTFds: Lowpass filter set to flp=%0.2f*sample_rate\n',default_flp);
  bandwidthMessage=1;
end
ds=updateBandwidth(ds,fhp,flp);

%  Update date/time fields of infods and res4
ds=updateDateTime(ds);

%  Create the .res4 file in the output dataset.
ds.res4=writeRes4([path,baseName,'.ds',delim,baseName,'.res4'],ds.res4,MAX_COILS);

if ds.res4.numcoef<0
  fprintf('\nwriteCTFds: writeRes4 returned ds.res4.numcoef=%d (<0??)\n\n',...
    ds.res4.numcoef);
  %  Kill the output dataset.
  rmdir([path,baseName,'.ds'],'s');
  ds=-1;
  return
end

%  Create .hist file
histfile=[path,baseName,'.ds',delim,baseName,'.hist'];
fid=fopen(histfile,'w');
fwrite(fid,ds.hist,'uint8');
fclose(fid);

% New .newds file
if isfield(ds,'newds')
  fid=fopen([path,baseName,'.ds',delim,baseName,'.newds'],'w');
  fwrite(fid,ds.newds,'uint8');
  fclose(fid);
end

%  New infods file.
if isfield(ds,'infods')
  writeCPersist([path,baseName,'.ds',delim,baseName,'.infods'],ds.infods);
end

% new hc file
if isfield(ds,'hc')
  ds.hc=writeHc([path,baseName,'.ds',delim,baseName,'.hc'],ds.hc,HLCdata(:,:,1));
end
clear HLCdata;

% MarkerFile.mrk
if checkMrk(ds)
  writeMarkerFile([path,baseName,'.ds',delim,'MarkerFile.mrk'],ds.mrk);
end

% .EEG
if isfield(ds,'EEG')
  writeEEG([path,baseName,'.ds',delim,baseName,'.EEG'],ds.EEG);
end

% .acq
if isfield(ds,'acq');
  % Check that ds.acq has the correct fields
  if isfield(ds.acq,'name') & isfield(ds.acq,'type') & isfield(ds.acq,'data')
    acqFilename=[path,baseName,'.ds',delim,baseName,'.acq'];
    writeCPersist(acqFilename,ds.acq);
  end
end

%  bad.segments file
if isfield(ds,'badSegments')
  writeBadSegments([path,baseName,'.ds',delim,'bad.segments'],ds.badSegments,...
    ds.res4.no_trials,ds.res4.no_samples/ds.res4.sample_rate);
end

% BadChannels
if isfield(ds,'BadChannels');
  if ischar(ds.BadChannels) & ~isempty(ds.BadChannels)
    fid=fopen([path,baseName,'.ds',delim,'BadChannels'],'w','ieee-be');
    for k=1:size(ds.BadChannels,1)
      fprintf(fid,'%s\n',deblank(ds.BadChannels(k,:)));
    end
    fclose(fid);
  end
end

%  ClassFile
if check_cls(ds)
  writeClassFile([path,baseName,'.ds',delim,'ClassFile.cls'],ds.TrialClass);
end

% VirtualChannels
if isfield(ds,'Virtual')
  writeVirtualChannels([path,baseName,'.ds',delim,'VirtualChannels'],ds.Virtual);
end

% processing.cfg
if isfield(ds,'processing');
  if ischar(ds.processing)
    fid=fopen([datasetname,delim,'processing.cfg'],'w','ieee-be');
    fwrite(fid,ds.processing,'uint8');
    fclose(fid);
  end
end

%  Update the data set path and name
ds.path=path;
ds.baseName=baseName;
return
% *************** End of function writeCTFds ********************************************
% **************************************************************************************


% **************************************************************************************
% *************** Function check_senres **************************
function [senres,status]=check_senres(senres,numChan);
% A user may have augmented the senres table, so check that all the fields have the
% correct size.    This will cause errors in programs that try to compute 
% sensor response using the geometry in the senres table.

%  Does "sanity" checks on the senres table.  If there are obviously incorrect entries, it
%  tries to fix them, and prints a message. If the senres table does not specify  coil 
%  positions, orientations or areas,  set them to zero, but give them the correct array
%  size.
newChannelType=4;  % Create the fake sensors as MEG magnetometers.  This way offsets can be removed.
status=-1;
%  Does the senres table have the correct no. of channels?
no_senres=length(senres);
if no_senres<numChan
  % Add channels.  Assume that they are gradiometers.
  ioGain=1;
  qGain=2^20;
  gain=0.3;
  properGain=1e15/(qGain*ioGain*gain);  % sensor gain in phi0/fT
  for kx=(no_senres+1):numChan
    senres(kx)=struct(...
      'sensorTypeIndex',newChannelType,'originalRunNum',0,'coilShape',0,...
      'properGain',properGain,'qGain',qGain,'ioGain',ioGain,...
      'ioOffset',0,'numCoils',1,...
      'grad_order_no',0,...
      'gain',gain,...
      'pos0',zeros(3,1),'ori0',zeros(3,1),...
      'area',0.1,'numturns',1,...
      'pos',[0 0 21]','ori',zeros(3,1));
  end
  fprintf(['\tAdded %d SQUID channels to senres table.  Nominal gain of each = ',...
      '%8.4f fT/step\n'],numChan-no_senres,gain);
  clear kx gain qGain ioGain properGain;
elseif no_senres>numChan
  %  Channels have been removed from the data set, but can't tell which elements of the
  %  senres table to remove.
  fprintf(['length(senres)=%d, but only %d channels of data.  writeCTFds can''t',...
      ' tell which channels to remove from senres table.\n'],no_senres,numChan);
  return
end
no_senres=length(senres);

%  Previous version of check_senres ensures that several fields had size [1 1].
%  Seems unnecessary, so removed it.

for k=1:no_senres
  % Check the fields that define pickup loops.  Force the field defining the pickup
  % loops to have the correct size.  It is not clear why the EEG channels and 
  % the type=13,28,29 HLC channels need numCoils=1, and area>0.
  if any(senres(k).sensorTypeIndex==[0:7])   % SQUID channels
    correct_numCoils=rem(senres(k).sensorTypeIndex,4)+1;
  elseif any(senres(k).sensorTypeIndex==[13 28 29])  % HLC channels
    correct_numCoils=1;
  elseif any(senres(k).sensorTypeIndex==[8 9])  % EEG channels
    correct_numCoils=1;
  else
    correct_numCoils=0;
  end
  if senres(k).numCoils~=correct_numCoils & any(senres(k).sensorTypeIndex==[0:7])
    fprintf('writeCTFds_test: senres(%d).sensorTypeIndex=%d but numCoils=%d??\n',...
      k,senres(k).sensorTypeIndex,senres(k).numCoils);
    fprintf('               Set numCoils=%d\n',correct_numCoils);
    senres(k).numCoils=correct_numCoils;
  end
  numCoils=senres(k).numCoils;
  if size(senres(k).pos0)~=[3 numCoils];pos0=zeros(3,numCoils);end
  if size(senres(k).ori0)~=[3 numCoils];ori0=zeros(3,numCoils);end
  if size(senres(k).pos)~=[3 numCoils];pos=zeros(3,numCoils);end
  if size(senres(k).ori)~=[3 numCoils];ori=zeros(3,numCoils);end
  if size(senres(k).numturns)~=[1 numCoils];numturns=zeros(1,numCoils);end
  if size(senres(k).area)~=[1 numCoils];area=zeros(3,numCoils);end
end
status=1;
return
%  ************* End of function check_senres*********************************************
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


% **************************************************************************************
%  ************* Function writeHc*********************************************
function hc=writeHc(hcFileName,hc,HLCdata);

%  Modified for v1.2 ds.hc structures.  ds.hc.names has the exact, complete names of the
%  head coils.  The coordinates relative to the subject may be ds.hc.head(MEG) OR
%  ds.hc.abdomen (fMEG).

%  Creates a .hc file in a CTF dataset
% Inputs: hcFileName : Complete name of .hc file including path, basename and .hc ext.
%         hc : structure with nasion, left, right head coil positions in dewar and
%              CTF head coordinates.
%         HLCdata : head coil locations at the start of 1st trial.  unit=m
%                   coordinates=dewar.  Used only if structure hc is empty and it is MEG
%                   data (not fMEG data).

% Output : hc : hc structure.  hc=struct([]) on failure.  hc is set to the coil positions
%                              in array HLCdata if hc=[] on entry.

%  Check inputs
if exist('hc')~=1
  hc=struct([]);
elseif ~isstruct(hc)
  hc=struct([]);
end
hcfields=fieldnames(hc);

if exist('HLCdata')~=1
  HLCdata=[];
elseif ~isnumeric(HLCdata) | size(HLCdata,1)~=3 | size(HLCdata,2)<3
  HLCdata=[];
end

%  Check hc Filename
if exist('hcFileName')~=1;hcFileName=char([]);end
if ~ischar(hcFileName) | isempty(hcFileName)
  fprintf('writeCTFds (writeHc): Called writeHc with bad file name.\n');
  hcFileName
  return
end

% Both hc and HLCdata bad?
if length(hcfields)~=4 & isempty(HLCdata)
  fprintf('writeCTFds (writeHc): Called writeHc with bad hc and bad HLCdata.\n');
  hc=struct([])
  return
elseif length(hcfields)~=4
  rstandard=8/sqrt(2)*[1 1 0;-1 1 0;1 -1 0]';
  rstandard(3,:)=-27;
  rdewar=100*HLCdata(:,1:3);
  
  %  Convert from dewar coordinates to CTF head coordinates.
  originCTF=0.5*(hc.dewar(:,2)+hc.dewar(:,3));
  % Unit vectors for the CTF coordinates
  uCTF(:,1)=hc.dewar(:,1)-originCTF;
  uCTF(:,3)=cross(uxCTF,hc.dewar(:,2)-hc.dewar(:,3));
  uCTF(:,2)=cross(uCTF(:,3),uCTF(:,1));
  uCTF=uCTF./(ones(3,1)*sqrt(sum(uCTF.^2,1)));
  rCTF=uCTF'*(rdewar-originCTF*ones(1,3))
  hc=struct('names',strvcat('nasion','left ear','right ear'),...
    'standard',rstandard,'dewar',rdewar,'head',rCTF);
  clear originCTF uCTF rstandard rdewar rCTF;
end

%  Character strings for generating the .hc text file
%  Should never have both 'head' and 'abdomen' fields.
labelword=strvcat('standard','measured','measured');
printField=[strmatch('standard',hcfields) strmatch('dewar',hcfields) ...
    strmatch('head',hcfields) strmatch('abdomen',hcfields)];

if ~strmatch('names',hcfields) | length(printField)~=3
  fprintf(['writeCTFds (writeHc): Structure hc does not have all of the required fields.\n',...
      '                           No .hc file will appear in the output dataset.\n']);
  hc;
  hc=struct([]);
  return
end

relative=strvcat('dewar','dewar',hcfields{printField(3)});
coilcoord=strvcat('standard','dewar',hcfields{printField(3)});

comp='xyz';
coilname=hc.names;

fid=fopen(hcFileName,'w','ieee-be');

for k=1:size(coilcoord,1)
  rcoil=getfield(hc,coilcoord(k,:));
  for coil=1:size(hc.names,1)
    clName=deblank(hc.names(coil,:));
    fwrite(fid,[labelword(k,:) ' ' clName ' coil position relative to ',...
        deblank(relative(k,:)) ' (cm):' char(10)],'uint8');
    for m=1:3
      fwrite(fid,[char(9) comp(m) ' = ' num2str(rcoil(m,coil),'%7.5f') char(10)],'uint8');
    end
  end
end
fclose(fid);
status=0;
return  
%  ************* End of function writeHc*********************************************
% **************************************************************************************


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Function checkMrk  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MarkerFileOK=checkMrk(ds);

%  Examines structure ds to see if a sensible MarkerFile can be created from ds.mrk.
%  Output: MarkerFileOK=1 : ds.mrk looks OK
%          MarkerFileOK=0 : ds.mrk cannot be a set of valid markers for these data.

MarkerFileOK=isfield(ds,'mrk');
if MarkerFileOK
  MarkerFileOK=~isempty(ds.mrk);
end
if MarkerFileOK
  %  Are the markers appropriate?
  minMarkerTrial=[];
  minMarkerTime=[];
  maxMarkerTrial=[];
  maxMarkerTime=[];
  for k=1:length(ds.mrk)
    maxMarkerTrial=max([maxMarkerTrial max(ds.mrk(k).trial)]);
    maxMarkerTime=max([maxMarkerTime max(ds.mrk(k).time)]);
    minMarkerTrial=min([minMarkerTrial min(ds.mrk(k).trial)]);
    minMarkerTime=min([minMarkerTime min(ds.mrk(k).time)]);
  end
  if isempty(maxMarkerTrial) | isempty(maxMarkerTime)
    MarkerFileOK=0;  % Do not create MarkerFile.mrk if all of the marker classes are empty.
  else
    MarkerFileOK=(maxMarkerTrial<=ds.res4.no_trials & minMarkerTrial>=1 & ...
      maxMarkerTime<=(ds.res4.no_samples/ds.res4.sample_rate) & ...
      minMarkerTime>=(-ds.res4.preTrigPts/ds.res4.sample_rate));
    if ~MarkerFileOK
      fprintf(['\nwriteCTFds (checkMrk): ds.mrk cannot possibly be a set of markers ',...
          'for array(data).\n']);
      fprintf(['          minMarkerTrial=%d (must be >=1)  ',...
          'maxMarkerTrial=%d (must be <=%d)\n'],...
        minMarkerTrial,maxMarkerTrial,ds.res4.no_trials);
      fprintf(['          minMarkerTime=%7.4f (must be >=%7.4f)  ',...
          'maxMarkerTrial=%7.4f  (must be <=%7.4f )\n'],...
        minMarkerTime,-ds.res4.preTrigPts/ds.res4.sample_rate,...
        maxMarkerTime,ds.res4.no_samples/ds.res4.sample_rate);
      fprintf('          MarkerFile.mrk will not be created.\n\n');
    end
  end
end
return
%%%%%%%%%%%%%% end of checkMrk %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function writeEEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeEEG(eegFileName,EEG);

%  Reads the EEG file of a dtaset and stores the infoemation in strucure array EEG.
%  EEG(k).chanName = channel name in the dataset (EEGmmm where mmm=channel number)
%  EEG(k).name = channel name given by the user (e.g. Fp4)
%  EEG(k).pos = electrode position in cm in CTF head coordinates
%  Check inputs

if exist('eegFileName')~=1;eegFileName=char([]);end
if ~ischar(eegFileName) | isempty(eegFileName)
  fprintf('writeCTFds (writeEEG): Called writeEEG with bad file name.\n');
  eegFileName
  EEG=struct([]);
end

if exist('EEG')~=1
  EEG=struct([]);
elseif ~isstruct(EEG)
  EEG=struct([]);
end

if isempty(EEG);return;end

fid=fopen(eegFileName,'w','ieee-be');
if fid<0
  fprintf('writeCTFds (writeEEG): Could not open file %s\n',eegFileName);
  return
end
nEEG=length(EEG);
for k=1:nEEG
  fprintf(fid,'%d\t%s\t%7.5f\t%7.5f\t%7.5f\n',EEG(k).chanNum,EEG(k).name,EEG(k).pos);
end
fclose(fid);
return
%%%%%%%%%  End of writeEEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Function writeBadSegments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeBadSegments(badSegmentsFile,badSegments,nTrial,tTrial);

%  Creates a bad.segements file in a CTF data set from the information in structure
%  badSegments which is created by read_badSegments, or by the user.
%  If structure badSegments is empty, then the file is not created.

%  badSegments structure:
%    badSegments.trial = List of trial numbers   
%    badSegments.StartTime = List of bad segment start times (relative to trial).
%    badSegments.EndTime   = List of bad segment end times.

%  Check badSegmentsFile
if exist('badSegmentsFile')~=1;badSegmentsFile=char([]);end
if isempty(badSegmentsFile) | ~ischar(badSegmentsFile)
  fprintf('writeCTFds(writeBadSegments):  Bad file name.\n');
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
  fprintf(['\nwriteCTFds (writeBadSegments): ',...
      'The fields of structure badSegments do not all have the same size.\n']);
  return
elseif any(badSegments.trial>nTrial) | any(badSegments.trial<1) | ...
    any(badSegments.EndTime>tTrial)
  fprintf(['\nwriteCTFds (writeBadSegments): ds.badSegments cannot possibly describe ',...
      'bad segments for these data.\n',...
      '\tmin(badSegments.trial)=%d max(badSegments.trial)=%d ',...
      'max(badSegments.EndTime)=%0.4f s\n\t\tDataset: nTrial=%d  tTrial=%0.4f s\n'],...
    min(badSegments.trial),max(badSegments.trial),max(badSegments.EndTime),...
    nTrial,tTrial);
  fprintf('\t\tbad.segments file will not be created.\n\n');
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
%%%%%%%%%%Function check_cls  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClassFileOK=check_cls(ds);

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
      fprintf(['\nwriteCTFds (check_cls): ds.TrialClass cannot possibly be a set of ',...
          'trial classes for array(data).\n          minClassTrial=%d (must be >=1)  ',...
          'maxClassTrial=%d (must be <=%d)\n'],...
        minClassTrial,maxClassTrial,ds.res4.no_trials);
      fprintf('          ClassFile.cls will not be created.\n');
    end
  end
end
return
%%%%%%%%%%%%%% end of check_cls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Function writeVirtualChannels    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeVirtualChannels(VirtualChannelsFile,Virtual);

%  Writes a VirtualChannels file using the information in structure Virtual.
%  S\tructure Virtual is prepared by read_VirtualChannels

%  Check VirtualChannelsFile
if exist('VirtualChannelsFile')~=1;VirtualChannelsFile=char([]);end
if isempty(VirtualChannelsFile) | ~ischar(VirtualChannelsFile)
  fprintf('write_VirtualChannelsFile:  Bad file name.\n');
  VirtualChannelsFile
  return
end

%  Check structure array Virtual
if exist('Virtual')~=1
  fprintf('writeVirtualChannels:  Must specify structure Virtual.\n');
  return
elseif isempty(Virtual) | ~isstruct(Virtual)
  return
elseif ~isfield(Virtual,'Name') | ~isfield(Virtual,'wt');
  return
elseif isempty(Virtual(1).Name)
  return
end

fid=fopen(VirtualChannelsFile,'w','ieee-be');
if fid<0
  fprintf('writeCTFds (writeVirtualChannels): Could not open file %s\n',VirtualChannelsFile);
  return
end
fprintf(fid,'//Virtual channel configuration\n\n');
for k=1:length(Virtual)
  fprintf(fid,'VirtualChannel\n{\n\tName:\t%s\n\tUnit:\t%s\n',...
    Virtual(k).Name,Virtual(k).Unit);
  for q=1:size(Virtual(k).wt)
    %  Floating format chosen to match VirtualChanels file creatd by 
    %  DataEditor (5.3.0-experimental-linux-20060918). 
    fprintf(fid,'\tRef:\t%s,%0.6g\n',deblank(Virtual(k).chan(q,:)),Virtual(k).wt(q));
  end
  fprintf(fid,'}\n\n');
end
fclose(fid);
return
%%%%%%%%%%%%%% end of writeVirtualChannels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Function updateDescriptors   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ds=updateDescriptors(ds,comment,creatorSoftware,HLCcntrl);

%  Makes sure that certain tags are in structure infods, and makes sure that information
%  is transfered to infods and newds from res4.  Updates several fields of res4 with
%  clinical use message (=comment) and name of creator program.

%  Inputs: ds : ds structure produced by readCTFds.
%          comment : Character string that is added to infods tags listed in addCommentTag
%              and res4 fields listed in addCommentField
%          creatorSoftware : Character string indicating that the data set was
%              created by writeCTFds.  Added to infods tags listed in addCreatorTag and
%              appName field of res4.
%          HLCcntrl: If ds.infods is missing or empty, HLCcntrl determines the
%             _DATA_HZ_MODE tage of infods.  if noit present or empty, HLCcntrl=0.
%  Creates a dummy infods structure if necessary using function make_dummy_infods.

%  Changes infods and newds to match information in 
%     ds.res4.nf_run_title
%     ds.res4.collect_descriptor
%     ds.res4.nf_subject_id
%     ds.res4.nf_operator

%  Adds comment (clinical use message) and creator software name to infods, newds and hist files.

if ~isfield(ds,'infods');ds.infods=[];end


%  infods tags that will display the comment.
addCommentTag=strvcat('_PATIENT_ID','_PATIENT_INSTITUTE','_PROCEDURE_COMMENTS',...
  '_DATASET_COMMENTS','_DATASET_STATUS','_DATASET_PROCSTEPTITLE');
%  res4 text fields that will display the comment
addCommentField=strvcat('appName','dataOrigin','dataDescription',...
  'nf_run_title','nf_subject_id','run_description');
%  res4 text string lengths (from document "CTF MEG File Formats", PN900-0088) 
addCommentLength=[256 256 256 256 32 -1]';  % -1 indicates variable length

%  infods tags that will display the creator software.
addCreatorTag=strvcat('_PROCEDURE_COMMENTS','_DATASET_STATUS',...
  '_DATASET_COLLECTIONSOFTWARE','_DATASET_CREATORSOFTWARE');
%  res4 text fields that will display the creator software.
addCreatorField=strvcat('appName','run_description');
addCreatorLength=[256 -1]';

%  List of infods tags that will be set to ds.res4.fields (not necessarily fields listed in
%  addCommentField or addCreatorField
addRes4Tag=strvcat('_PATIENT_ID','_DATASET_PROCSTEPTITLE',...
  '_DATASET_PROCSTEPDESCRIPTION','_DATASET_OPERATORNAME');

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Add the Tags to infods.  In most cases they will already be there, but just to be sure.
addTag=strvcat(addRes4Tag,addCommentTag,addCreatorTag);
tagClass=strvcat('_PATIENT','_PROCEDURE','_DATASET');
% List of all the tag names as a character array
tagName=getArrayField(ds.infods,'name')';
if length(ds.infods)==1;tagName=tagName';end
tagType=getArrayField(ds.infods,'type');
tagPtr=1;  % If a class is missing, inject the class starting at tagPtr+1.

for k=1:size(tagClass,1)
  addIndex=strmatch(deblank(tagClass(k,:)),addTag)';
  % Are there any tags to be added in this class?
  if isempty(addIndex);continue;end
  % List of infods tags in the tagClass, but excluding the CPerist type (which marks the
  % start of a class.
  infodsIndex=strmatch(deblank(tagClass(k,:)),tagName);
  if isempty(infodsIndex)  % Create a new class of tags.
    if strcmp(deblank(tagName(tagPtr+1,:)),'EndOfParameters') & k>1;
      tagPtr=tagPtr+1;
    end
    nTag=length(ds.infods);
    tagName((tagPtr+4):(nTag+3),:)=tagName((tagPtr+1):nTag,:);
    ds.infods((tagPtr+4):(nTag+3))=ds.infods((tagPtr+1):nTag);
    ds.infods(tagPtr+1)=struct('name',[deblank(tagClass(k,:)),'_INFO'],'type',2,'data',[]);
    ds.infods(tagPtr+2)=struct('name','WS1_','type',0,'data',[]);
    ds.infods(tagPtr+3)=struct('name','EndOfParameters','type',-1,'data',[]);
    tagName=strvcat(tagName(1:tagPtr,:),...
      strvcat([deblank(tagClass(k,:)),'_INFO'],'WS1_','EndOfParameters'),...
      tagName(tagPtr+1:nTag,:));
    nTag=nTag+3;
    tagPtr=tagPtr+2;
  else
    if ds.infods(max(infodsIndex)).type==2
      tagPtr=max(infodsIndex)+1; % Class consists of no tags beyond the CPersist definition
    else
      tagPtr=max(infodsIndex);   % Class contains at least one real tag.
    end
    clear infodsIndex;
  end
  for q=addIndex
    if isempty(strmatch(deblank(addTag(q,:)),tagName))
      tagName=strvcat(tagName(1:tagPtr,:),deblank(addTag(q,:)),tagName(tagPtr+1:nTag,:));
      ds.infods((tagPtr+2):(nTag+1))=ds.infods((tagPtr+1):nTag);
      ds.infods(tagPtr+1)=struct('name',deblank(addTag(q,:)),'type',10,'data',char([]));
      tagPtr=tagPtr+1;
      nTag=nTag+1;
    end
  end
end
clear k q nTag tagPtr infodsIndex addIndex;
%  All of the tags in addTag have been added to ds.infods.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('creatorSoftware')~=1;
  creatorSoftware=char([]);
elseif ~ischar(creatorSoftware)
  creatorSoftware=char([]);
else
  creatorSoftware=deblank(creatorSoftware);
end

if exist('comment')~=1;
  comment=char([]);
elseif ~ischar(comment)
  comment=char([]);
else
  comment=deblank(comment);
end

tagName=getArrayField(ds.infods,'name')';

%  Update the res4 fields : add creator message.
%  Don't check field length, but truncate later.
for q=1:size(addCreatorField,1);
  strng=deblank(getfield(ds.res4,deblank(addCreatorField(q,:))));
  if isempty(strfind(strng,creatorSoftware))
    newStrng=creatorSoftware;
    if ~isempty(strng);newStrng=[strng '  ' newStrng];end
    ds.res4=setfield(ds.res4,deblank(addCreatorField(q,:)),newStrng);
  end
end
clear q strng nChar strng0 ns;

%  Update the res4 fields : add comment (= clinical use message)
%  Don't check field length, but truncate later.
for q=1:size(addCommentField,1);
  strng=deblank(getfield(ds.res4,deblank(addCommentField(q,:))));
  if isempty(strfind(strng,comment))
    newStrng=comment;
    if ~isempty(strng);newStrng=[strng '  ' newStrng];end
    ds.res4=setfield(ds.res4,deblank(addCommentField(q,:)),newStrng);
  end
end
clear strng newStrng q;

% Update res4.run_description
ds.res4.run_description=[ds.res4.run_description,char(0)];
ds.res4.rdlen=length(ds.res4.run_description);

%  Update infods entries from .res4 fields
ds.infods(strmatch('_PATIENT_ID',tagName)).data=deblank(ds.res4.nf_subject_id);
ds.infods(strmatch('_DATASET_PROCSTEPTITLE',tagName)).data=deblank(ds.res4.nf_run_title);
ds.infods(strmatch('_DATASET_PROCSTEPDESCRIPTION',tagName)).data=...
  deblank(ds.res4.nf_collect_descriptor);
ds.infods(strmatch('_DATASET_OPERATORNAME',tagName)).data=deblank(ds.res4.nf_operator);

%  Truncate the .res4. fields.  Leave room for a final char(0).
for q=1:size(addCreatorField,1);
  strng=deblank(getfield(ds.res4,deblank(addCreatorField(q,:))));
  if length(strng)>addCreatorLength(q)-1 & addCreatorLength(q)>0
    ds.res4=setfield(ds.res4,deblank(addCreatorField(q,:)),...
      strng(length(strng)+[-addCreatorLength(q)+2:0]));
  end
end
for q=1:size(addCommentField,1);
  strng=deblank(getfield(ds.res4,deblank(addCommentField(q,:))));
  if length(strng)>addCommentLength(q)-1 & addCommentLength(q)>0
    ds.res4=setfield(ds.res4,deblank(addCommentField(q,:)),...
      strng(length(strng)+[-addCommentLength(q)+2:0]));
  end
end
clear q strng;

% Add creator software to infods tags.  Have already cheked that the tags are there.
for q=1:size(addCreatorTag,1)
  if isempty(strmatch(deblank(addCreatorTag(q,:)),addRes4Tag))
    k=strmatch(deblank(addCreatorTag(q,:)),tagName);
    if length(k)==1
      if isempty(strfind(ds.infods(k).data,creatorSoftware))
        newStrng=creatorSoftware;
        if ~isempty(deblank(ds.infods(k).data));
          newStrng=[deblank(ds.infods(k).data) '  ' newStrng];
        end
        ds.infods(k).data=newStrng;
      end
    else
      fprintf('writeCTFds: Tag %s appears %d times in ds.infods ??\n',...
        deblank(addCreatorTag(q,:)),length(k));
    end
  end
end
clear q k;

% Add comment (clinical use statement) to ds.infods
for q=1:size(addCommentTag,1)
  if isempty(strmatch(deblank(addCommentTag(q,:)),addRes4Tag))
    k=strmatch(deblank(addCommentTag(q,:)),tagName);
    if length(k)==1
      if isempty(strfind(ds.infods(k).data,comment))
        newStrng=comment;
        if ~isempty(deblank(ds.infods(k).data));
          newStrng=[deblank(ds.infods(k).data) '  ' newStrng];
        end
        ds.infods(k).data=newStrng;
      end
    else
      fprintf(['writeCTFds (updateDescriptors): Tag %s appears %d times in ',...
          'ds.infods ??\n'],deblank(addCommentTag(q,:)),length(k));
    end
  end
end
clear q k;

%  Add the creator and comment information to the newds file
if isfield(ds,'newds');
  nChar=length(ds.newds);
  %  Keyword positions
  aNpos=max(strfind(ds.newds,[char(9) 'appName:' char(9)]));
  if isempty(aNpos)
    fprintf(['writeCTFds (update_decsriptors): Keyword ''appName'' ',...
        'does not appear in ds.newds.\n',...
        '         set ds.newds=[].\n']);
    ds.newds=char([]);
  else
    eol=max([11 min(strfind(ds.newds(aNpos+1:length(ds.newds)),char(10)))]);
    ds.newds=[ds.newds(1:(aNpos+eol-1)) '  ' creatorSoftware '  ' comment ...
        ds.newds((aNpos+eol):length(ds.newds))];
  end
  clear nChar aNpos eol;
end

% Add clinical use message to hist file.
if isfield(ds,'hist') & ~isempty(comment)
  ds.hist=[ds.hist char(10) comment char(10)];
end

return
%%%%%%%%%%%%%% End of updateDescriptors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Function updateHLC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [infods,status]=updateHLC(infods,HLCdata);
%  Adds notes and  continuous-head-localization tags to .infods file

%  Inputs : infods : Structure of tags to be written to output infods file.
%           HLCdata : 0 or [] or not defined - There are no HLC channels
%                     (head-coil position) in the data.  Just pass the existing tags to
%                     the output dataset.
%           HLCdata : real array : size(HLCdata)=[samplesPerTrial 3 Ncoil nTrial]
%                     Calculate the head motion tags.

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

%  Check HLCdata array.  HLCdata=[] indicates no continuous head localization.
if ~exist('HLCdata')==1
  HLCdata=[];
elseif ~isempty(HLCdata)
  [samplesPerTrial nDim Ncoil nTrial]=size(HLCdata);
  if nDim~=3;
    fprintf('writeCTFds (updateHLC): size(HLCdata)=[');fprintf(' %d',size(HLCdata));
    fprintf(']\n');
    fprintf('             nDim=%d??   Setting HLCdata=[].\n',nDim);
    HLCdata=[];
  end
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
    HLCoffset=squeeze(HLCdata(1,:,:,1));
    for q=1:Ncoil
      for k=1:3
        HLCdata(:,k,q,:)=HLCdata(:,k,q,:)-HLCoffset(k,q);  
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
    TagValue(3)=maxHeadMotion;             % _DATASET_MAXHEADMOTION
    %  Subtract 1 from the trial so trial numbering starts at 0.
    TagValue(4)=maxHeadMotionTrial-1;      % _DATASET_MAXHEADMOTIONTRIAL
    TextTagValue=char(zeros(4,1)); % Placeholder since tags1:4 are numerical
    
    % _MAXHEADMOTIONCOIL value
    TextTagValue=strvcat(TextTagValue,sprintf('%d',maxHeadMotionCoil));
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  Function updateBandwidth  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ds=updateBandwidth(ds,fhp,flp);
%  Updates the bandwidth fields of infods and newds.
%  3 possible dources of bandwidth:
%     1. fhp,flp inputs   (default)
%     2. newds fiels
%     3. res4.filter
%     4. infods field
%  Almost always, fhp,flp will have the default values, but it could be defined in the
%  fields of ds, so check.  Result: lots of code to achieve a simple result.

%  read fhp,flp from newds (if it exists)
if isfield(ds,'newds')
  BWpos=max(strfind(ds.newds,[char(9) 'bandwidth:' char(9)])); %  Keyword position
  if isempty(BWpos)
    fprintf(['writeCTFds (updateBandwidth): Keyword ''bandwidth:''',...
        'does not appear in ds.newds.\n',...
        '         Field newds removed from structure ds.\n']);
    ds=rmfield(ds,'newds');
  else
    %  Get the bandwidth parameters from ds.newds.
    eol=max([13 min(strfind(ds.newds(BWpos+1:length(ds.newds)),char(10)))]);
    buff=sscanf(ds.newds((BWpos+13):(BWpos+eol)),'%f%c%f');
    fhp=max(fhp,buff(1));
    flp=min(flp,buff(3));
  end
  clear eol buff BWpos;
end

%  Get fhp, flp from res4.
if ~isempty(ds.res4.filters)
  for kq=1:ds.res4.num_filters
    freq=ds.res4.filters(kq).freq;
    if ds.res4.filters(kq).fType==1;
      flp=min(flp,freq);
    elseif ds.res4.filters(kq).fType==2;
      fhp=max(fhp,freq);
    end
  end
  clear kq freq;
end

%  Get fhp, flp from ds.infods
if isfield(ds,'infods')
  name=getArrayField(ds.infods,'name')';
  TagNo=strmatch('_DATASET_LOWERBANDWIDTH',name,'exact');
  if ~isempty(TagNo);fhp=max(fhp,ds.infods(TagNo).data);end
  TagNo=strmatch('_DATASET_UPPERBANDWIDTH',name,'exact');
  if ~isempty(TagNo);flp=min(flp,ds.infods(TagNo).data);end
end

%  Now have fhp,flp.  Update hist.  Add all the res4 filters.
if ~isempty(ds.res4.filters)
  if ~isfield(ds,'hist');hist=char([]);end
  ds.hist=[ds.hist 'Filters specified in res4 file :' char(10)];
  for kq=1:ds.res4.num_filters
    freq=ds.res4.filters(kq).freq;
    if ds.res4.filters(kq).fType==1;
      ds.hist=[ds.hist ...
          num2str(kq,'%8d') '  Lowpass ' num2str(freq,'%6.2f') ' Hz' char(10)];
    elseif ds.res4.filters(kq).fType==2;
      ds.hist=[ds.hist ...
          num2str(kq,'%8d') '  Highpass ' num2str(freq,'%6.2f') ' Hz' char(10)];
    elseif ds.res4.filters(kq).fType==3;
      ds.hist=[ds.hist num2str(kq,'%8d') '  Notch ' num2str(freq,'%6.2f')  ...
          ' Hz    width=' num2str(ds.res4.filters(kq).Param,'%6.2f') ' Hz' char(10)];
    else
      ds.hist=[ds.hist '   fType=',num2str(ds.res4.filters(kq).fType,'%d') char(10)];
    end
  end
  ds.hist=[ds.hist ,...
      'Bandwidth: ',num2str(fhp,'%0.3g'),' - ',num2str(flp,'%0.3g'),' Hz',char(10)];
  clear kq freq;
end

% Update newds bandwidth
if isfield(ds,'newds')
  BWpos=max(strfind(ds.newds,[char(9) 'bandwidth:' char(9)]));
  eol=max([12 min(strfind(ds.newds(BWpos+1:length(ds.newds)),char(10)))]);
  ds.newds=[ds.newds(1:BWpos+11) num2str(fhp,'%0.6f') ', ' num2str(flp,'%0.6f') ...
      ds.newds((BWpos+eol):length(ds.newds))];
  clear BWpos eol;
end

%  Update infods bandwidth tags.  Lots of coding becasue it's possible that infods does
%  not already have bandwidth tags.
name=getArrayField(ds.infods,'name')';
DATASET_tags=strmatch('_DATASET_',name);
if ~isempty(DATASET_tags)
  addTag=strvcat('_DATASET_LOWERBANDWIDTH','_DATASET_UPPERBANDWIDTH');
  TagValue=[fhp flp];
  % Add or update tags.
  for q=1:size(addTag,1)
    nTag=length(ds.infods);
    TagNo=strmatch(deblank(addTag(q,:)),name,'exact')';
    if isempty(TagNo)  % Insert a tag at the end of the _DATASET tags.
      TagNo=max(DATASET_tags)+1;
      ds.infods((TagNo+1):(nTag+1))=ds.infods(TagNo:nTag);
      ds.infods(TagNo)=struct('name',deblank(addTag(q,:)),'type',4,'data',TagValue(q));
      name=strvcat(name(1:TagNo-1,:),deblank(addTag(q,:)),name(TagNo:nTag,:));
      DATASET_tags=[DATASET_tags TagNo];
    else
      ds.infods(TagNo).data=TagValue(q);  % Tag exists.  Just update the value.
    end
  end  % End loop over head position and head motion tags.
  clear q nTag DATASET_tags;
else
  fprintf(['writeCTFds (updateBandwidth):  Found no _DATASET tags in infods.\n'...
      '\t\t\tDid not add bandwidth tags.\n']);
end

return
%%%%%%%%%%%%%% End of updateBandwidth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  Function updateDateTime  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ds=updateDateTime(ds);
%  Update the dat/time fields of res4 and and infods

datetime=floor(clock);
if isfield(ds,'infods');
  name=getArrayField(ds.infods,'name')';
  DATASET_tags=strmatch('_DATASET_',name);
  if ~isempty(DATASET_tags)
    addTag=strvcat('_DATASET_CREATORDATETIME','_DATASET_LASTMODIFIEDDATETIME');
    tagData=sprintf('%d%2.2d%2.2d%2.2d%2.2d%2.2d',datetime);
    tagType=10;
    for q=1:size(addTag,1)
      nTag=length(ds.infods);
      TagNo=strmatch(deblank(addTag(q,:)),name,'exact')';
      if isempty(TagNo)  % Insert a tag at the end of the _DATASET tags.
        TagNo=max(DATASET_tags)+1;
        ds.infods((TagNo+1):(nTag+1))=ds.infods(TagNo:nTag);
        name=strvcat(name(1:TagNo-1,:),deblank(addTag(q,:)),name(TagNo:nTag,:));
        DATASET_tags=[DATASET_tags TagNo];
      end
      ds.infods(TagNo)=struct('name',deblank(addTag(q,:)),'type',tagType,'data',tagData);
    end
    clear nTag TagNo name DATASET_tags addTag tagData tagType;
  else
    fprintf('writeCTFds (updateDateTime): Cannot find any _DATASET tags in ds.infods.\n');
  end  % End loop over head position and head motion tags.
else
  fprintf('writeCTFds (updateDateTime): Cannot find ds.infods.\n');
end
ds.res4.data_date=sprintf('%4d/%2.2d/%2.2d',datetime(1:3));
ds.res4.data_time=sprintf('%2.2d:%2.2d:%2.2d',datetime(4:6));

return
%%%%%%%%%%%%%% End of updateDateTime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
