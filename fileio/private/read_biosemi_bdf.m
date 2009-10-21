function dat = read_biosemi_bdf(filename, hdr, begsample, endsample, chanindx);

% READ_BIOSEMI_BDF reads specified samples from a BDF continous datafile
% It neglects all trial boundaries as if the data was acquired in
% non-continous mode.
%
% Use as
%   [hdr] = read_biosemi_bdf(filename);
% where
%    filename        name of the datafile, including the .bdf extension
% This returns a header structure with the following elements
%   hdr.Fs           sampling frequency
%   hdr.nChans       number of channels
%   hdr.nSamples     number of samples per trial
%   hdr.nSamplesPre  number of pre-trigger samples in each trial
%   hdr.nTrials      number of trials
%   hdr.label        cell-array with labels of each channel
%   hdr.orig         detailled EDF header information
%
% Or use as
%   [dat] = read_biosemi_bdf(filename, hdr, begsample, endsample, chanindx);
% where
%    filename        name of the datafile, including the .bdf extension
%    hdr             header structure, see above
%    begsample       index of the first sample to read
%    endsample       index of the last sample to read
%    chanindx	     index of channels to read (optional, default is all)
% This returns a Nchans X Nsamples data matrix

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: read_biosemi_bdf.m,v $
% Revision 1.2  2009/10/12 12:14:23  roboos
% fixed small typo in readLowLevel fnuction, thanks to Philip
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.7  2008/10/06 08:42:17  roboos
% added missing "end" in subfunction
%
% Revision 1.6  2008/10/03 12:20:57  roboos
% use matlab reading instead of mex file in case file >2GB, thanks to Philip
%
% Revision 1.5  2008/09/30 07:47:04  roboos
% replaced all occurences of setstr() with char(), because setstr is deprecated by Matlab
%
% Revision 1.4  2007/10/01 13:44:21  roboos
% changed a detail in the calibration (in case of one channel the output would remain sparse)
%
% Revision 1.3  2007/09/13 09:52:20  roboos
% removed section of code regarding idx1/2/3 which was slow but unneeded
% avoid confusion between EDF and hdr.orig
% implemented much faster reading for a single channel (efficient when reading status channel)
%
% Revision 1.2  2007/09/12 12:52:09  roboos
% calibrate the data
%
% Revision 1.1  2006/02/01 08:18:48  roboos
% new implementation, based on some EEGLAB code and a mex file for the binary part of the data
%

if nargin==1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the header, this code is from EEGLAB's openbdf
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FILENAME = filename;

  % defines Seperator for Subdirectories
  SLASH='/';
  BSLASH=char(92);

  cname=computer;
  if cname(1:2)=='PC' SLASH=BSLASH; end;

  fid=fopen(FILENAME,'r','ieee-le');
  if fid<0
    fprintf(2,['Error LOADEDF: File ' FILENAME ' not found\n']);
    return;
  end;

  EDF.FILE.FID=fid;
  EDF.FILE.OPEN = 1;
  EDF.FileName = FILENAME;

  PPos=min([max(find(FILENAME=='.')) length(FILENAME)+1]);
  SPos=max([0 find((FILENAME=='/') | (FILENAME==BSLASH))]);
  EDF.FILE.Ext = FILENAME(PPos+1:length(FILENAME));
  EDF.FILE.Name = FILENAME(SPos+1:PPos-1);
  if SPos==0
    EDF.FILE.Path = pwd;
  else
    EDF.FILE.Path = FILENAME(1:SPos-1);
  end;
  EDF.FileName = [EDF.FILE.Path SLASH EDF.FILE.Name '.' EDF.FILE.Ext];

  H1=char(fread(EDF.FILE.FID,256,'char')');     %
  EDF.VERSION=H1(1:8);                          % 8 Byte  Versionsnummer
  %if 0 fprintf(2,'LOADEDF: WARNING  Version EDF Format %i',ver); end;
  EDF.PID = deblank(H1(9:88));                  % 80 Byte local patient identification
  EDF.RID = deblank(H1(89:168));                % 80 Byte local recording identification
  %EDF.H.StartDate = H1(169:176);               % 8 Byte
  %EDF.H.StartTime = H1(177:184);               % 8 Byte
  EDF.T0=[str2num(H1(168+[7 8])) str2num(H1(168+[4 5])) str2num(H1(168+[1 2])) str2num(H1(168+[9 10])) str2num(H1(168+[12 13])) str2num(H1(168+[15 16])) ];

  % Y2K compatibility until year 2090
  if EDF.VERSION(1)=='0'
    if EDF.T0(1) < 91
      EDF.T0(1)=2000+EDF.T0(1);
    else
      EDF.T0(1)=1900+EDF.T0(1);
    end;
  else ;
    % in a future version, this is hopefully not needed
  end;

  EDF.HeadLen = str2num(H1(185:192));  % 8 Byte  Length of Header
  % reserved = H1(193:236);	           % 44 Byte
  EDF.NRec = str2num(H1(237:244));     % 8 Byte  # of data records
  EDF.Dur = str2num(H1(245:252));      % 8 Byte  # duration of data record in sec
  EDF.NS = str2num(H1(253:256));       % 8 Byte  # of signals

  EDF.Label = char(fread(EDF.FILE.FID,[16,EDF.NS],'char')');
  EDF.Transducer = char(fread(EDF.FILE.FID,[80,EDF.NS],'char')');
  EDF.PhysDim = char(fread(EDF.FILE.FID,[8,EDF.NS],'char')');

  EDF.PhysMin= str2num(char(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));
  EDF.PhysMax= str2num(char(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));
  EDF.DigMin = str2num(char(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));
  EDF.DigMax = str2num(char(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));

  % check validity of DigMin and DigMax
  if (length(EDF.DigMin) ~= EDF.NS)
    fprintf(2,'Warning OPENEDF: Failing Digital Minimum\n');
    EDF.DigMin = -(2^15)*ones(EDF.NS,1);
  end
  if (length(EDF.DigMax) ~= EDF.NS)
    fprintf(2,'Warning OPENEDF: Failing Digital Maximum\n');
    EDF.DigMax = (2^15-1)*ones(EDF.NS,1);
  end
  if (any(EDF.DigMin >= EDF.DigMax))
    fprintf(2,'Warning OPENEDF: Digital Minimum larger than Maximum\n');
  end
  % check validity of PhysMin and PhysMax
  if (length(EDF.PhysMin) ~= EDF.NS)
    fprintf(2,'Warning OPENEDF: Failing Physical Minimum\n');
    EDF.PhysMin = EDF.DigMin;
  end
  if (length(EDF.PhysMax) ~= EDF.NS)
    fprintf(2,'Warning OPENEDF: Failing Physical Maximum\n');
    EDF.PhysMax = EDF.DigMax;
  end
  if (any(EDF.PhysMin >= EDF.PhysMax))
    fprintf(2,'Warning OPENEDF: Physical Minimum larger than Maximum\n');
    EDF.PhysMin = EDF.DigMin;
    EDF.PhysMax = EDF.DigMax;
  end
  EDF.PreFilt= char(fread(EDF.FILE.FID,[80,EDF.NS],'char')');	%
  tmp = fread(EDF.FILE.FID,[8,EDF.NS],'char')';	%	samples per data record
  EDF.SPR = str2num(char(tmp));               %	samples per data record

  fseek(EDF.FILE.FID,32*EDF.NS,0);

  EDF.Cal = (EDF.PhysMax-EDF.PhysMin)./(EDF.DigMax-EDF.DigMin);
  EDF.Off = EDF.PhysMin - EDF.Cal .* EDF.DigMin;
  tmp = find(EDF.Cal < 0);
  EDF.Cal(tmp) = ones(size(tmp));
  EDF.Off(tmp) = zeros(size(tmp));

  EDF.Calib=[EDF.Off';(diag(EDF.Cal))];
  %EDF.Calib=sparse(diag([1; EDF.Cal]));
  %EDF.Calib(1,2:EDF.NS+1)=EDF.Off';

  EDF.SampleRate = EDF.SPR / EDF.Dur;

  EDF.FILE.POS = ftell(EDF.FILE.FID);
  if EDF.NRec == -1                            % unknown record size, determine correct NRec
    fseek(EDF.FILE.FID, 0, 'eof');
    endpos = ftell(EDF.FILE.FID);
    EDF.NRec = floor((endpos - EDF.FILE.POS) / (sum(EDF.SPR) * 2));
    fseek(EDF.FILE.FID, EDF.FILE.POS, 'bof');
    H1(237:244)=sprintf('%-8i',EDF.NRec);      % write number of records
  end;

  EDF.Chan_Select=(EDF.SPR==max(EDF.SPR));
  for k=1:EDF.NS
    if EDF.Chan_Select(k)
      EDF.ChanTyp(k)='N';
    else
      EDF.ChanTyp(k)=' ';
    end;
    if findstr(upper(EDF.Label(k,:)),'ECG')
      EDF.ChanTyp(k)='C';
    elseif findstr(upper(EDF.Label(k,:)),'EKG')
      EDF.ChanTyp(k)='C';
    elseif findstr(upper(EDF.Label(k,:)),'EEG')
      EDF.ChanTyp(k)='E';
    elseif findstr(upper(EDF.Label(k,:)),'EOG')
      EDF.ChanTyp(k)='O';
    elseif findstr(upper(EDF.Label(k,:)),'EMG')
      EDF.ChanTyp(k)='M';
    end;
  end;

  EDF.AS.spb = sum(EDF.SPR);	% Samples per Block

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % convert the header to Fieldtrip-style
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if any(EDF.SampleRate~=EDF.SampleRate(1))
    error('channels with different sampling rate not supported');
  end
  hdr.Fs          = EDF.SampleRate(1);
  hdr.nChans      = EDF.NS;
  hdr.label       = cellstr(EDF.Label);
  % it is continuous data, therefore append all records in one trial
  hdr.nTrials     = 1;
  hdr.nSamples    = EDF.NRec * EDF.Dur * EDF.SampleRate(1);
  hdr.nSamplesPre = 0;
  hdr.orig        = EDF;

  % return the header
  dat = hdr;

else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % retrieve the original header
  EDF = hdr.orig;

  % determine the trial containing the begin and end sample
  epochlength = EDF.Dur * EDF.SampleRate(1);
  begepoch    = floor((begsample-1)/epochlength) + 1;
  endepoch    = floor((endsample-1)/epochlength) + 1;
  nepochs     = endepoch - begepoch + 1;
  nchans      = EDF.NS;

  if nargin<5
    chanindx = 1:nchans;
  end

  % allocate memory to hold the data
  dat = zeros(length(chanindx),nepochs*epochlength);

  % read and concatenate all required data epochs
  for i=begepoch:endepoch
    offset = EDF.HeadLen + (i-1)*epochlength*nchans*3;
    if length(chanindx)==1
      % this is more efficient if only one channel has to be read, e.g. the status channel
      offset = offset + (chanindx-1)*epochlength*3;
      buf = readLowLevel(filename, offset, epochlength); % see below in subfunction
      dat(:,((i-begepoch)*epochlength+1):((i-begepoch+1)*epochlength)) = buf;
    else
      % read the data from all channels and then select the desired channels
      buf = readLowLevel(filename, offset, epochlength*nchans); % see below in subfunction
      buf = reshape(buf, epochlength, nchans);
      dat(:,((i-begepoch)*epochlength+1):((i-begepoch+1)*epochlength)) = buf(:,chanindx)';
    end
  end

  % select the desired samples
  begsample = begsample - (begepoch-1)*epochlength;  % correct for the number of bytes that were skipped
  endsample = endsample - (begepoch-1)*epochlength;  % correct for the number of bytes that were skipped
  dat = dat(:, begsample:endsample);

  % Calibrate the data
  if length(chanindx)>1
    % using a sparse matrix speeds up the multiplication
    calib = sparse(diag(EDF.Cal(chanindx,:)));
    dat   = calib * dat;
  else
    % in case of one channel the calibration would result in a sparse array
    calib = diag(EDF.Cal(chanindx,:));
    dat   = calib * dat;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading the 24 bit values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buf = readLowLevel(filename, offset, numwords);
if offset < 2*1024^3
  % use the external mex file, only works for <2GB
  buf = read_24bit(filename, offset, numwords);
  % this would be the only difference between the bdf and edf implementation
  % buf = read_16bit(filename, offset, numwords);
else
  % use plain matlab, thanks to Philip van der Broek
  fp = fopen(filename,'r','ieee-le');
  status = fseek(fp, offset, 'bof');
  if status
    error(['failed seeking ' filename]);
  end
  [buf,num] = fread(fp,numwords,'bit24=>double');
  fclose(fp);
  if (num<numwords)
    error(['failed opening ' filename]);
    return
  end
end

