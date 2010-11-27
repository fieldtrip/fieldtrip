function data=getCTFdata(ds,dataList,chanList,unit,prec);

%  getCTFdata.m   Reads specified trials from .meg4 files in CTF-format data set.

% ************************************************************************
%
%   This program is provided to users of CTF MEG systems as a courtesy only.
%   It's operation has not been verified for clinical use.
%   Please do not redistribute it without permission from VSM MedTech Ltd.
%
% ************************************************************************

%  Version: 1.1  13 April 2007
%  Author : Harold Wilson

%  Reads trials(s) from CTF data sets and converts to physical units.
%  Modified to handle cases where the amount of data is so large (>2^31-8 bytes) 
%  that additional meg4 files ( .1_meg4,.2_meg4,...) are required.  

%  getCTFdata does not apply the fact that Acq creates multiple .meg4 files that are the
%  largest number of trials that fit into 2GB.  Rather, getCTFdata simply checks all of the
%  MEG4 files and makes a list of which trials are in which file.

%  Inputs : ds : Structure containing all of the res4 information on a file.
%                getCTFdata uses dtructure ds to get the meg4 filename , trial-size
%                information (ds.res4.no_samples, ds.res4.no_channels, ds.res4.no_trials)
%                and the gain information contained in the sensor resource table 
%                ds.res4.senres.

%     dataList : Tells getCTFdata which data to read.  Trial numbering starts at 1.
%           (1) dataList<=0 or [].  Read the complete data set.
%           (2) size(dataList)=[1 N].  Read the N trials listed in row vector dataList.
%           (3) size(dataList)=[3 1].  Read points dataList(1):dataList(2) 
%                                      from trial dataList(3).

%          unit :  Determines the unit of the signals:
%                   'ft' or 'fT' : Convert to fT (MEG), uV (EEG) and V (ADC) (default)
%                   't'  or 'T'  : Convert to  T (MEG),  V (EEG) and V (ADC)
%                   'phi0' : Convert to phi0 (MEG), uV (EEG) and V (ADC)
%                   'int': Read plain integers from *.meg4-file

%          chanList: List of channels to read.  Default : read all channels.
%             Options : (1) chanList<=0 or [].  Read all the channels
%                       (2) List of channel numbers.  Channel numbers are
%                           referred to the data set, NOT just the SQUID channels.
%                           Fail if any(chanList>ds.res4.no_channels)
%                       (3) A list of channel names.  chanList(n,:) is the name of
%                           the nth channel to be returned in array data.  getCTFdata
%                           compares the characters preceding the first blank or '-'
%                           character to the rows of ds.res4.chanNames.
%                           Fail if any channel name dows not appear in ds.res4.chanNames.
%                             

%          prec :  Precision of output.  Default='double'.
%                  User may select 'single' to reduce memory demands.

%  Output : data : Array of data.  Precision = double or single.
%                  size(data)=[nSample nChan nTrial] where
%                     nChan=ds.res4.no_channels or length(chanList) or size(chanList,1)
%                  If an error is detected, data is set to [];

%  Call to getCTFdata without arguments prints command format.
if nargout==0 & nargin==0
  fprintf(['\ngetCTFData: Version 1.1  13 April 2007   ',...
      'Reads a CTF dataset into array data.\n',...
      '\tCall : data=getCTFdata(ds,dataList,chanList,unit,prec);\n'...
      '\t\t   - dataList: <=0, [] or missing.  Read the complete data set.\n',...
      '\t\t                size(dataList)=[1 N].  Read the N trials listed.\n',...
      '\t\t                dataList=[pt1 pt2 Tr]''.  Read points pt1:pt2',...
      ' of trial Tr.\n']);
  return
end

%  Allowed 8-byte headers for meg4 files.
allowedMeg4Headers=strvcat(['MEG41CP',char(0)],['MEG42CP',char(0)]);

data=[];    % If an error causes in an early return.
delim=filesep;

if exist('ds')~=1
  fprintf('\ngetCTFdata: No input arguments.\n\n');
  return
elseif isempty(ds)
  fprintf('\ngetCTFdata: ds is empty.\n\n');
  return
elseif ~isstruct(ds)
  fprintf('\ngetCTFdata: Input class(ds)=%s.  Must be a structure.\n\n',class(ds));
  return
end 

%  Create the list of trials and range of points to read.
%  Set dataList=-1 to read all trials.
if exist('dataList')~=1;
  dataList=-1;
elseif isempty(dataList)
  dataList=-1;
elseif ~isnumeric(dataList)
  fprintf('getCTFdata: Input dataList has class=%s.  Must be numeric.\n',class(dataList));
  return
elseif ((any(dataList<=0) | any(dataList)>ds.res4.no_trials) & length(dataList)>1) |...
    sum(size(dataList)>1)>1
  fprintf('\ngetCTFdata: Size(dataList)=[');fprintf(' %d',size(dataList));fprintf(']\n');
  fprintf('\t\tmin(dataList)=%d  max(dataList)=%d   nTrial=%d\n',...
    min(reshape(dataList,1,prod(size(dataList)))),...
    max(reshape(dataList,1,prod(size(dataList)))),ds.res4.no_trials);
  fprintf(['\t\tdata_list must have size=[1 N] or [3 1].   ',...
      'No elements<=0 when dataList is not scalar.\n\n']);
  return
elseif isequal(size(dataList),[1 1])
  if dataList<=0;dataList=-1;end
end

if isequal(dataList,-1)
  trials=1:ds.res4.no_trials;    % read all trials
  point=[1 ds.res4.no_samples];  % read all points from each trial
elseif size(dataList,1)==1   % dataList is a list of trials
  if max(dataList)>ds.res4.no_trials | min(dataList)<=0
    fprintf(['\ngetCTFdata: min(dataList)=%d max(dataList)=%d.  ',...
        'ds.res4.no_trials=%d ??\n\n'],min(dataList),max(dataList),ds.res4.no_trials);
    return
  end
  trials=dataList;    % read listed trials
  point=[1 ds.res4.no_samples]; % read all points from each trial
elseif isequal(size(dataList),[3 1]) % dataList specifies points from a single trial
  if dataList(1)<=0 | dataList(1)>dataList(2) | dataList(2)>ds.res4.no_samples | ...
      dataList(3)<=0 | dataList(3)>ds.res4.no_trials
    fprintf('getCTFdata: dataList error: points=%d:%d  trial=%d\n',dataList);
    fprintf('            points/trial=%d   No. of trials=%d\n',...
      ds.res4.no_samples,ds.res4.no_trials);
    return
  end
  point=dataList(1:2)';
  trials=dataList(3);      % read only one trial
else
  fprintf('\ngetCTFdata: size(dataList)=[');fprintf(' %d',size(dataList));
  fprintf('] ??\n\n');
  return
end
clear dataList;

%  Check the specified signal unit.  Convert unit to lower case.
if exist('unit')~=1
  unit='ft';       % default
elseif isempty(unit)
  unit='ft';       % default
elseif ~ischar(unit);
  fprintf('\ngetCTFdata: Input unit has class=%s.  Must be char.\n\n',class(unit));
  return
else
  unit=lower(unit);
  if ~strcmp(unit,'int') & ~strcmp(unit,'ft') & ~strcmp(unit,'t') & ~strcmp(unit,'phi0')
    fprintf(['\ngetCTFdata : unit=%s    Not a valid option.  Must be ',...
        '''fT'', ''T'', ''phi0'' or ''int''\n\n'],unit);
    return
  end
end
  
%  Make a list of channels to read.
if exist('chanList')~=1;chanList=-1;end
if isempty(chanList);chanList=-1;end
if isnumeric(chanList)    % Channels listed by number
  chanList=reshape(chanList,1,prod(size(chanList)));
  if isequal(size(chanList),[1 1]) & chanList<=0;  % chanList<=0
    chanList=[1:ds.res4.no_channels];
  elseif max(chanList)>ds.res4.no_channels | min(chanList)<1  % error
    fprintf(['\ngetCTFdata: chanList specifies channels %d,%d but data set ',...
        'has only %d channels.\n\n'],min(chanList),max(chanList),ds.res4.no_channels);
    return
  end
elseif ischar(chanList)  % Channels listed by name
  name_list=chanList;
  nChan=size(name_list,1);
  chanList=zeros(1,nChan);
  for k=1:nChan
    kx=strmatch(strtok(name_list(k,:),'- '),ds.res4.chanNames);
    if isempty(kx)
      fprintf(['\ngetCTFdata: Channel %s does not appear in the list of channel ',...
          'names  (See ds.res4.chanNames)\n\n'],deblank(name_list(k,:)));
      return
    end
    chanList(k)=kx;
  end
  clear name_list kx k;
else
  fprintf('getCTFdata: Cannot interpret chanList.\n');
  return
end

%  Sort the channels into increasing number, and keep track of the original order.
%  C_list = Sorted list of channels read into array data.  The sorted list is used so
%           the meg4 file is read in increaing channel number order.
%  chan_index : C_list(p)=chanList(chan_index(p)).
%  data(:,k,:) is from dataset channel chanList(k).

[C_list chan_index]=sort(chanList);
nChan=length(chanList);

%  data array precision
if exist('prec')~=1;
  prec='double';
elseif isempty(prec)
  prec='double';
elseif ~ischar(prec)
  fprintf('\ngetCTFdata: Input prec has class=%s\n\n',class(prec));
  return
elseif ~strcmp(prec,'single') & ~strcmp(prec,'double')
  fprintf('\ngetCTFdata : prec=%s     Must be ''single'' or ''double''.\n\n',prec);
  return
end

% Data are stored in the meg4 file as 4-byte words.
trial_size=ds.res4.no_samples*ds.res4.no_channels;
bytes_per_trial=4*trial_size;
bytes_per_channel=4*ds.res4.no_samples;

nSavepts=point(2)-point(1)+1;
ntot=length(trials)*nChan*nSavepts;  % Total points to read.

%  Make a list of which trials are in which .meg4 file.  In almost all cases, only the .meg4
%  file is present.
qFile=-1;
while 1
  qFile=qFile+1;
  if qFile==0
    meg4Ext='.meg4';
    meg4Trial=0;
  else
    meg4Ext=['.',int2str(qFile),'_meg4'];
  end
  D=dir([ds.path,ds.baseName,'.ds',delim,ds.baseName,meg4Ext]);
  if isempty(D);break;end
  % No of trials in this file.
  nFileTrial=round((D.bytes-8)/bytes_per_trial);
  if nFileTrial*bytes_per_trial~=D.bytes-8
    fprintf('\ngetCTFdata: File %s%s does not contain a whole number of trials.\n',...
      ds.baseName,meg4Ext);
    fprintf(['            File size (bytes)=%0.9d    header=8 bytes   trial size ',...
        '(bytes)=%d\n\n'],D.bytes,bytes_per_trial);
    return
  end
  meg4Trial=[meg4Trial round((D.bytes-8)/bytes_per_trial)];
end

% meg4Trial(2)=last trial in file .meg4, meg4Trial(3)=last trial in file .1_meg4,...
meg4Trial=cumsum(meg4Trial);
nMEG4=length(meg4Trial)-1;

% Open .meg4 with big-endian option for compatibility with the data-acquisition computer.
meg4Fid=fopen([ds.path,ds.baseName,'.ds',delim,ds.baseName,'.meg4'],'r','ieee-be');
openFile=1;

% Check that the header of the .meg4 file is OK.  Don't bother to check the other
% meg4 files (.1_meg4,...)
meg4Header=char(fread(meg4Fid,8,'uint8')');
if isempty(strmatch(meg4Header,allowedMeg4Headers,'exact'))
  fclose(meg4Fid);
  fprintf(['\ngetCTFdata: meg4 file header = %s.  ',...
      'A valid meg4 file, must have one of'],meg4Header);
  for k=1:size(allowedMeg4Headers,1)
    fprintf(' %s',allowedMeg4Headers(k,:));
  end
  fprintf('\n\n');
  return
end

% Reserve memory for the data.
if prec=='single'   % Funny coding to minimize the maximum memory required.
  data=single([]);
  while length(data)<ntot
    data=[data single(zeros(1,min(round(1e7),ntot-length(data))))];
  end
else
  data=zeros(1,ntot);
end

data=reshape(data,nSavepts,nChan,length(trials));

for kt=1:length(trials)
  %  Check if this trial is in the open meg4 file.  If not open the correct meg4 file.
  if trials(kt)<=meg4Trial(openFile) | trials(kt)>meg4Trial(openFile+1)
    fclose(meg4Fid);
    openFile=min(find(meg4Trial>=trials(kt)))-1;
    if isempty(openFile)
      fprintf(['\ngetCTFdata: Trying to read trial %d, but the last trial in ',...
          'the dataset is %d.\n',...
          '            Dataset = %s\n\n'],trials(kt),max(meg4Trial),ds.baseName);
      data=[];
      return
    elseif openFile==1
      meg4Ext='.meg4';
    else
      meg4Ext=['.',int2str(openFile-1),'_meg4'];
    end
    meg4Fid=fopen([ds.path,ds.baseName,'.ds',delim,ds.baseName,meg4Ext],'r','ieee-be');
  end
  firstTrial=meg4Trial(openFile)+1;
  for kch=1:nChan
    % Move to the start of trial kt, channel C_list(kch), point(1)
    fseek(meg4Fid,8+(trials(kt)-firstTrial)*bytes_per_trial+...
      (C_list(kch)-1)*bytes_per_channel+(point(1)-1)*4,'bof');
    % Read and save nSavepts of data.
    if prec=='double'
      data(:,chan_index(kch),kt)=fread(meg4Fid,nSavepts,'int32');
    else
      data(:,chan_index(kch),kt)=single(fread(meg4Fid,nSavepts,'int32'));
    end
  end
end
fclose(meg4Fid);
clear kt kch meg4Fid;
clear openFile meg4Ext firstTrial;

if strcmp(unit,'int');return;end

% Convert to specified signal unit.
for kx=1:nChan
  sens=chanList(kx);
  SQUIDtype=any(ds.res4.senres(sens).sensorTypeIndex==[0:7]);
  EEGtype=any(ds.res4.senres(sens).sensorTypeIndex==[8 9]);
  if SQUIDtype & strcmp(unit,'phi0')
    invgain=ds.res4.senres(sens).ioGain*ds.res4.senres(sens).qGain;
  else
    invgain=ds.res4.senres(sens).ioGain*...
      ds.res4.senres(sens).properGain*ds.res4.senres(sens).qGain;
  end
  if strcmp(unit,'ft') & SQUIDtype
    pGain=1e15;  % Convert SQUIDs to fT
  elseif (strcmp(unit,'phi0') | strcmp(unit,'ft')) & EEGtype
    pGain=1e6;   % Convert from EEG to uV
  else
    pGain=1;
  end
  %  Trap case with invgain==0 which may be used to mark bad channels.
  if invgain~=0
    gainx=pGain/invgain;
  else
    gainx=0;
  end
  if abs(gainx-1)>1e-7
    for kt=1:length(trials)
      if prec=='single'
        data(:,kx,kt)=single(gainx*double(data(:,kx,kt)));
      else
        data(:,kx,kt)=gainx*data(:,kx,kt);
      end
    end
  end
end
clear kx kt sens SQUIDtype EEGtype pGain invGain gainx;

return
