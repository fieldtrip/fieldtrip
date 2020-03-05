function [dat] = read_edf(filename, hdr, begsample, endsample, chanindx)

% READ_EDF reads specified samples from an EDF datafile. It neglects all trial or
% data block boundaries as if the data was acquired in non-continuous mode.
%
% Note that since FieldTrip only accommodates a single sampling rate in a given
% dataset, whereas EDF allows specification of a sampling rate for each channel.  If
% there are heterogenous sampling rates then this function will automatically choose
% a subset.  If the last such channel is different from the rest, the assumption will
% be made that it is the annotation channel and the rest will be selected.  If that
% is not the case, then the largest subset of channels with a consistent sampling
% rate will be chosen.  To avoid this automatic selection process, the user may
% specify their own choice of channels using chanindx.  In this case, the automatic
% selection will only occur if the user selected channels still have heterogenous
% sampling rates.  In this case the automatic selection will occur amongst the user
% specified channels.  While reading the header the resulting channel selection
% decision will be stored in hdr.orig.chansel and the contents of this field will
% override chanindx during data reading.
%
% Use as
%   [hdr] = read_edf(filename)
% where
%    filename        name of the datafile, including the .edf extension
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
%   [dat] = read_edf(filename, hdr, begsample, endsample, chanindx)
% where
%    filename        name of the datafile, including the .edf extension
%    hdr             header structure, see above
%    begsample       index of the first sample to read
%    endsample       index of the last sample to read
%    chanindx        index of channels to read (optional, default is all)
% This returns a Nchans X Nsamples data matrix
%
% Or use as
%   [evt] = read_edf(filename, hdr)
% where
%    filename        name of the datafile, including the .edf extension
%    hdr             header structure, see above
% This returns an Nsamples data vector of just the annotation channel

% Copyright (C) 2006, Robert Oostenveld and others
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

switch nargin
  case 1
    chanindx=[];
  case 2
    chanindx=[];
  case 3
    chanindx=begsample;
  case 4
end

needhdr = (nargin==1)||(nargin==3);
needevt = (nargin==2);
needdat = (nargin==5);

if needhdr
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the header, this code is from EEGLAB's openbdf
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FILENAME = filename;

  % defines Seperator for Subdirectories
  SLASH='/';
  BSLASH=char(92);

  cname=computer;
  if cname(1:2)=='PC' SLASH=BSLASH; end

  fid=fopen_or_error(FILENAME,'r','ieee-le');

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
  end
  EDF.FileName = [EDF.FILE.Path SLASH EDF.FILE.Name '.' EDF.FILE.Ext];

  H1=char(fread(EDF.FILE.FID,256,'char')');
  EDF.VERSION=H1(1:8);                          % 8 Byte  Versionsnummer
  %if 0 fprintf(2,'LOADEDF: WARNING  Version EDF Format %i',ver); end
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
    end
  else
    % in a future version, this is hopefully not needed
  end

  EDF.HeadLen = str2num(H1(185:192));  % 8 Byte  Length of Header
  % reserved = H1(193:236);            % 44 Byte
  EDF.NRec = str2num(H1(237:244));     % 8 Byte  # of data records
  EDF.Dur  = str2num(H1(245:252));     % 8 Byte  # duration of data record in sec
  EDF.NS   = str2num(H1(253:256));     % 8 Byte  # of signals

  EDF.Label      = char(fread(EDF.FILE.FID,[16,EDF.NS],'char')');
  EDF.Transducer = char(fread(EDF.FILE.FID,[80,EDF.NS],'char')');
  EDF.PhysDim    = char(fread(EDF.FILE.FID,[ 8,EDF.NS],'char')');

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
    fprintf(2,'Warning OPENEDF: Failing Physical Minimum, taking Digital Minimum instead\n');
    EDF.PhysMin = EDF.DigMin;
  end
  if (length(EDF.PhysMax) ~= EDF.NS)
    fprintf(2,'Warning OPENEDF: Failing Physical Maximum, taking Digital Maximum instead\n');
    EDF.PhysMax = EDF.DigMax;
  end
  idx_PhysMin_ge_PhysMax = EDF.PhysMin >= EDF.PhysMax;
  if (any(idx_PhysMin_ge_PhysMax))
    tmplabel = cellfun(@(x) [x ' '], cellstr(EDF.Label(idx_PhysMin_ge_PhysMax,:)),'UniformOutput',false)';
    fprintf(2,['Warning OPENEDF: Physical Minimum larger than Maximum.\nPLEASE recheck if the scaling and polarity in the following channels are still correct if used:\n' tmplabel{:} '\n']);
    %EDF.PhysMin = EDF.DigMin;
    %EDF.PhysMax = EDF.DigMax;
  end
  EDF.PreFilt= char(fread(EDF.FILE.FID,[80,EDF.NS],'char')');
  EDF.SPR = str2num(char(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));  % samples per data record

  fseek(EDF.FILE.FID,32*EDF.NS,0);

  EDF.Cal = (EDF.PhysMax-EDF.PhysMin)./(EDF.DigMax-EDF.DigMin);
  EDF.Cal(isnan(EDF.Cal)) = 0;
  EDF.Off = EDF.PhysMin - EDF.Cal .* EDF.DigMin;
  
  %tmp = find(EDF.Cal < 0);
  %EDF.Cal(tmp) = ones(size(tmp));
  %EDF.Off(tmp) = zeros(size(tmp));

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
  end

  EDF.Chan_Select=(EDF.SPR==max(EDF.SPR));
  for k=1:EDF.NS
    if EDF.Chan_Select(k)
      EDF.ChanTyp(k)='N';
    else
      EDF.ChanTyp(k)=' ';
    end
    if contains(upper(EDF.Label(k,:)),'ECG')
      EDF.ChanTyp(k)='C';
    elseif contains(upper(EDF.Label(k,:)),'EKG')
      EDF.ChanTyp(k)='C';
    elseif contains(upper(EDF.Label(k,:)),'EEG')
      EDF.ChanTyp(k)='E';
    elseif contains(upper(EDF.Label(k,:)),'EOG')
      EDF.ChanTyp(k)='O';
    elseif contains(upper(EDF.Label(k,:)),'EMG')
      EDF.ChanTyp(k)='M';
    end
  end

  if isempty(chanindx)
    chanindx=1:EDF.NS;
  end

  EDF.AS.spb = sum(EDF.SPR);    % Samples per Block
  bi=[0;cumsum(EDF.SPR)];

  idx=[];idx2=[];
  for k=1:EDF.NS
    idx2=[idx2, (k-1)*max(EDF.SPR)+(1:EDF.SPR(k))];
  end
  maxspr=max(EDF.SPR);
  idx3=zeros(EDF.NS*maxspr,1);
  for k=1:EDF.NS, idx3(maxspr*(k-1)+(1:maxspr))=bi(k)+ceil((1:maxspr)'/maxspr*EDF.SPR(k));end

  %EDF.AS.bi=bi;
  EDF.AS.IDX2=idx2;
  %EDF.AS.IDX3=idx3;

  % close the file
  fclose(EDF.FILE.FID);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % convert the header to Fieldtrip-style
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if all(EDF.SampleRate(chanindx)==EDF.SampleRate(chanindx(1)))
    chansel=chanindx;
    hdr.Fs           = EDF.SampleRate(chanindx(1));
    hdr.nChans       = length(chansel);
    hdr.label        = cellstr(EDF.Label(chansel,:));
    % it is continuous data, therefore append all records in one trial
    hdr.nSamples     = EDF.NRec * EDF.SPR(chansel(1));
    hdr.nSamplesPre  = 0;
    hdr.nTrials      = 1;
    hdr.chanunit     = cellstr(EDF.PhysDim(chansel,:));
    hdr.chantype     = repmat({'unknown'}, size(hdr.chanunit));  % start with unknown
    hdr.chantype(strcmp(hdr.chanunit, 'uV')) = {'eeg'};          % it might also be EOG, ECG, EMG, etc
    hdr.chantype(strcmp(hdr.chanunit, 'Boolean')) = {'trigger'};
    hdr.orig         = EDF;
    % this will be used on subsequent reading of data
    if length(chansel) ~= EDF.NS
      hdr.orig.chansel = chansel;
    else
      hdr.orig.chansel = 1:hdr.nChans;
    end
    hdr.orig.annotation = find(strcmp(cellstr(hdr.orig.Label), 'EDF Annotations'));

  elseif all(EDF.SampleRate(1:end-1)==EDF.SampleRate(1))
    % only the last channel has a deviant sampling frequency
    % this is the case for EGI recorded datasets that have been converted
    % to EDF+, in which case the annotation channel is the last
    chansel = find(EDF.SampleRate==EDF.SampleRate(1));
    % continue with the subset of channels that has a consistent sampling frequency
    hdr.Fs           = EDF.SampleRate(chansel(1));
    hdr.nChans       = length(chansel);
    ft_warning('Skipping "%s" as continuous data channel because of inconsistent sampling frequency (%g Hz)', deblank(EDF.Label(end,:)), EDF.SampleRate(end));
    hdr.label        = cellstr(EDF.Label(chansel,:));
    % it is continuous data, therefore append all records in one trial
    hdr.nSamples     = EDF.NRec * EDF.SPR(chansel(1));
    hdr.nSamplesPre  = 0;
    hdr.nTrials      = 1;
    hdr.orig         = EDF;
    % this will be used on subsequent reading of data
    hdr.orig.chansel    = chansel;
    hdr.orig.annotation = find(strcmp(cellstr(hdr.orig.Label), 'EDF Annotations'));

  else
    % select the sampling rate that results in the most channels
    [a, b, c] = unique(EDF.SampleRate);
    chancount = nan(size(a));
    for i=1:length(a)
      chancount(i) = sum(c==i);
    end
    [dum, indx] = max(chancount);
    chansel = find(EDF.SampleRate == a(indx));

    % continue with the subset of channels that has a consistent sampling frequency
    hdr.Fs           = EDF.SampleRate(chansel(1));
    hdr.nChans       = length(chansel);
    hdr.label        = cellstr(EDF.Label);
    hdr.label        = hdr.label(chansel);
    % it is continuous data, therefore append all records in one trial
    hdr.nSamples     = EDF.NRec * EDF.SPR(chansel(1));
    hdr.nSamplesPre  = 0;
    hdr.nTrials      = 1;
    hdr.orig         = EDF;
    % this will be used on subsequent reading of data
    hdr.orig.chansel    = chansel;
    hdr.orig.annotation = find(strcmp(cellstr(hdr.orig.Label), 'EDF Annotations'));

    ft_warning('channels with different sampling rate not supported, selecting subset of %d channels at %f Hz', length(hdr.label), hdr.Fs);
  end

  % return the header
  dat = hdr;

elseif needdat || needevt
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % retrieve the original header
  EDF = hdr.orig;

  % There can be an optional chansel field containing a list of predefined channels.
  % These channels are in that case also the only ones represented in the FieldTrip
  % header, which means that the other channels are simply not visible to the naive
  % user. This field can be present because the user specified an explicit channel
  % selection in FT_READ_HEADER or because the read_edf function had to automatically
  % choose a subset to cope with heterogenous sampling rates or even both.  In any
  % case, at this point in the file reading process the contents of the chansel field
  % has the proper specification for channel selection, taking into account both the
  % user channel selection as well as any correction that might have been made due to
  % heterogenous sampling rates.

  if     ~isempty(chanindx) && ~isfield(EDF, 'chansel')
    % a subset of channels should been selected from the full list of channels in the file
    chanindx = chanindx; % keep as it is
    useChanindx = true;
  elseif ~isempty(chanindx) &&  isfield(EDF, 'chansel')
    % a subset of channels should been selected from the predefined list
    chanindx = EDF.chansel(chanindx);
    useChanindx = true;
  elseif  isempty(chanindx) &&  isfield(EDF, 'chansel')
    % all channels from the predefined list should be selected
    chanindx = EDF.chansel(chanindx);
    useChanindx = true;
  elseif  isempty(chanindx) && ~isfield(EDF, 'chansel')
    %  simply select all channels that are present in the file
    chanindx = 1:EDF.NS;
    useChanindx = false;
  end

  if needevt
    % read the annotation channel, not the data channels
    chanindx = EDF.annotation;
    begsample = 1;
    endsample = EDF.SPR(EDF.annotation)*EDF.NRec;
  end

  if useChanindx
    epochlength = EDF.SPR(chanindx(1));   % in samples for the selected channel
    blocksize   = sum(EDF.SPR);           % in samples for all channels
    chanoffset  = EDF.SPR;
    chanoffset  = round(cumsum([0; chanoffset(1:end-1)]));
    nchans      = length(chanindx);       % get the selection from the subset of channels
  else
    epochlength = EDF.SPR(1);             % in samples for a single channel
    blocksize   = sum(EDF.SPR);           % in samples for all channels
    nchans      = EDF.NS;                 % use all channels
  end

  % determine the trial containing the begin and end sample
  begepoch    = floor((begsample-1)/epochlength) + 1;
  endepoch    = floor((endsample-1)/epochlength) + 1;
  nepochs     = endepoch - begepoch + 1;

  % allocate memory to hold the data
  dat = zeros(length(chanindx),nepochs*epochlength);

  % read and concatenate all required data epochs
  for i=begepoch:endepoch
    if useChanindx
      % only a subset of channels with consistent sampling frequency is read
      offset = EDF.HeadLen + (i-1)*blocksize*2; % in bytes
      % read the complete data block
      buf = readLowLevel(filename, offset, blocksize); % see below in subfunction
      for j=1:length(chanindx)
        % cut out the part that corresponds with a single channel
        dat(j,((i-begepoch)*epochlength+1):((i-begepoch+1)*epochlength)) = buf((1:epochlength) + chanoffset(chanindx(j)));
      end

    elseif length(chanindx)==1
      % this is more efficient if only one channel has to be read, e.g. the status channel
      offset = EDF.HeadLen + (i-1)*blocksize*2; % in bytes
      offset = offset + (chanindx-1)*epochlength*2;
      % read the data for a single channel
      buf = readLowLevel(filename, offset, epochlength); % see below in subfunction
      dat(:,((i-begepoch)*epochlength+1):((i-begepoch+1)*epochlength)) = buf;

    else
      % read the data from all channels, subsequently select the desired channels
      offset = EDF.HeadLen + (i-1)*blocksize*2; % in bytes
      % read the complete data block
      buf = readLowLevel(filename, offset, blocksize); % see below in subfunction
      buf = reshape(buf, epochlength, nchans);
      dat(:,((i-begepoch)*epochlength+1):((i-begepoch+1)*epochlength)) = buf(:,chanindx)';
    end
  end

  % select the desired samples
  begsample = begsample - (begepoch-1)*epochlength;  % correct for the number of bytes that were skipped
  endsample = endsample - (begepoch-1)*epochlength;  % correct for the number of bytes that were skipped
  dat = dat(:, begsample:endsample);

  % Calibrate the data
  if useChanindx
    calib = diag(EDF.Cal(chanindx));
  end
  if length(chanindx)>1
    % using a sparse matrix speeds up the multiplication
    dat = sparse(calib) * dat;
  else
    % in case of one channel the sparse multiplication would result in a sparse array
    dat = calib * dat;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading the 16 bit values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buf = readLowLevel(filename, offset, numwords)
is_below_2GB = offset < 2*1024^2;
read_16bit_success = true;
if is_below_2GB
  % use the external mex file, only works for <2GB
  try
  buf = read_16bit(filename, offset, numwords);
  catch e
      read_16bit_success = false;
  end
end
if ~is_below_2GB || ~read_16bit_success
  % use plain matlab, thanks to Philip van der Broek
  fp = fopen(filename,'r','ieee-le');
  status = fseek(fp, offset, 'bof');
  if status
    ft_error(['failed seeking ' filename]);
  end
  [buf,num] = fread(fp,numwords,'bit16=>double');
  fclose(fp);
  if (num<numwords)
    ft_error(['failed reading ' filename]);
    return
  end
end
