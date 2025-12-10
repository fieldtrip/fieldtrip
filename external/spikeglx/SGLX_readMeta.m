% =============================================================
% Simple class of MATLAB functions deomonstrating
% how to read and manipulate SpikeGLX meta and binary files.
% When this file is included in the MATLAB path, the individual methods
% are called with SGLX_readMeta.<function name>
%
% The most important method in the class is ReadMeta().
% The functions for reading binary data are not optimized for speed,
% and are included as a useful starting point.
%
classdef SGLX_readMeta

  methods(Static)

    % =========================================================
    % Parse ini file returning a structure whose field names
    % are the metadata left-hand-side tags, and whose right-
    % hand-side values are MATLAB strings. We remove any
    % leading '~' characters from tags because MATLAB uses
    % '~' as an operator.
    %
    % If you're unfamiliar with structures, the benefit
    % is that after calling this function you can refer
    % to metafile items by name. For example:
    %
    %   meta.fileCreateTime  // file create date and time
    %   meta.nSavedChans     // channels per timepoint
    %
    % All of the values are MATLAB strings, but you can
    % obtain a numeric value using str2double(meta.nSavedChans).
    % More complicated parsing of values is demonstrated in the
    % utility functions below.
    %
    function [meta] = ReadMeta(binName, path)

      % Create the matching metafile name
      [~,name,~] = fileparts(binName);
      metaName = strcat(name, '.meta');

      % Parse ini file into cell entries C{1}{i} = C{2}{i}
      fid = fopen(fullfile(path, metaName), 'r');
      % -------------------------------------------------------------
      %    Need 'BufSize' adjustment for MATLAB earlier than 2014
      %    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
      C = textscan(fid, '%[^=] = %[^\r\n]');
      % -------------------------------------------------------------
      fclose(fid);

      % New empty struct
      meta = struct();

      % Convert each cell entry into a struct entry
      for i = 1:length(C{1})
        tag = C{1}{i};
        if tag(1) == '~'
          % remake tag excluding first character
          tag = sprintf('%s', tag(2:end));
        end
        meta.(tag) = C{2}{i};
      end
    end % ReadMeta


    % =========================================================
    % Read nSamp timepoints from the binary file, starting
    % at timepoint offset samp0. The returned array has
    % dimensions [nChan,nSamp]. Note that nSamp returned
    % is the lesser of: {nSamp, timepoints available}.
    %
    % IMPORTANT: samp0 and nSamp must be integers.
    %
    function dataArray = ReadBin(samp0, nSamp, meta, binName, path)

      nChan = str2double(meta.nSavedChans);

      nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
      samp0 = max(samp0, 0);
      nSamp = min(nSamp, nFileSamp - samp0);

      sizeA = [nChan, nSamp];

      fid = fopen(fullfile(path, binName), 'rb');
      fseek(fid, samp0 * 2 * nChan, 'bof');
      dataArray = fread(fid, sizeA, 'int16=>double');
      fclose(fid);
    end % ReadBin


    % =========================================================
    % Return an array [lines X timepoints] of uint8 values for
    % a specified set of digital lines.
    %
    % - dwReq is the one-based index into the saved file of the
    %    16-bit word that contains the digital lines of interest.
    % - dLineList is a zero-based list of one or more lines/bits
    %    to scan from word dwReq.
    %
    function digArray = ExtractDigital(dataArray, meta, dwReq, dLineList)
      % Get channel index of requested digital word dwReq
      if strcmp(meta.typeThis, 'imec')
        [AP, LF, SY] = SGLX_readMeta.ChannelCountsIM(meta);
        if SY == 0
          fprintf('No imec sync channel saved\n');
          digArray = [];
          return;
        else
          digCh = AP + LF + dwReq;
        end
      else
        [MN,MA,XA,DW] = SGLX_readMeta.ChannelCountsNI(meta);
        if dwReq > DW
          fprintf('Maximum digital word in file = %d\n', DW);
          digArray = [];
          return;
        else
          digCh = MN + MA + XA + dwReq;
        end
      end
      [~,nSamp] = size(dataArray);
      digArray = zeros(numel(dLineList), nSamp, 'uint8');
      for i = 1:numel(dLineList)
        digArray(i,:) = bitget(dataArray(digCh,:), dLineList(i)+1, 'int16');
      end
    end % ExtractDigital


    % =========================================================
    % Return sample rate as double.
    %
    function srate = SampRate(meta)
      if strcmp(meta.typeThis, 'imec')
        srate = str2double(meta.imSampRate);
      else
        srate = str2double(meta.niSampRate);
      end
    end % SampRate


    % =========================================================
    % Return a multiplicative factor for converting 16-bit
    % file data to voltage. This does not take gain into
    % account. The full conversion with gain is:
    %
    %   dataVolts = dataInt * fI2V / gain.
    %
    % Note that each channel may have its own gain.
    %
    function fI2V = Int2Volts(meta)
      if strcmp(meta.typeThis, 'imec')
        if isfield(meta,'imMaxInt')
          maxInt = str2num(meta.imMaxInt);
        else
          maxInt = 512;
        end
        fI2V = str2double(meta.imAiRangeMax) / maxInt;
      else
        fI2V = str2double(meta.niAiRangeMax) / 32768;
      end
    end % Int2Volts


    % =========================================================
    % Return array of original channel IDs. As an example,
    % suppose we want the imec gain for the ith channel stored
    % in the binary data. A gain array can be obtained using
    % ChanGainsIM() but we need an original channel index to
    % do the look-up. Because you can selectively save channels
    % the ith channel in the file isn't necessarily the ith
    % acquired channel, so use this function to convert from
    % ith stored to original index.
    %
    % Note: In SpikeGLX channels are 0-based, but MATLAB uses
    % 1-based indexing, so we add 1 to the original IDs here.
    %
    function chans = OriginalChans(meta)
      if strcmp(meta.snsSaveChanSubset, 'all')
        chans = (1:str2double(meta.nSavedChans));
      else
        chans = str2num(meta.snsSaveChanSubset);
        chans = chans + 1;
      end
    end % OriginalChans


    % =========================================================
    % Return counts of each imec channel type that compose
    % the timepoints stored in binary file.
    %
    function [AP,LF,SY] = ChannelCountsIM(meta)
      M = str2num(meta.snsApLfSy);
      AP = M(1);
      LF = M(2);
      SY = M(3);
    end % ChannelCountsIM

    % =========================================================
    % Return counts of each nidq channel type that compose
    % the timepoints stored in binary file.
    %
    function [MN,MA,XA,DW] = ChannelCountsNI(meta)
      M = str2num(meta.snsMnMaXaDw);
      MN = M(1);
      MA = M(2);
      XA = M(3);
      DW = M(4);
    end % ChannelCountsNI


    % =========================================================
    % Return gain for ith channel stored in the nidq file.
    %
    % ichan is a saved channel index, rather than an original
    % (acquired) index.
    %
    function gain = ChanGainNI(ichan, savedMN, savedMA, meta)
      if ichan <= savedMN
        gain = str2double(meta.niMNGain);
      elseif ichan <= savedMN + savedMA
        gain = str2double(meta.niMAGain);
      else
        gain = 1;
      end
    end % ChanGainNI


    % =========================================================
    % Return gain arrays for imec channels.
    %
    % Index into these with original (acquired) channel IDs.
    %
    function [APgain,LFgain] = ChanGainsIM(meta)
      % list of probe types with NP 1.0 imro format
      np1_imro = [0,1020,1030,1200,1100,1120,1121,1122,1123,1300];
      % number of channels acquired
      acqCountList = str2num(meta.acqApLfSy);

      APgain = zeros(acqCountList(1));     % default type = float64
      LFgain = zeros(acqCountList(2));     % empty array for 2.0

      if isfield(meta,'imDatPrb_type')
        probeType = str2double(meta.imDatPrb_type);
      else
        probeType = 0;
      end

      if ismember(probeType, np1_imro)
        % imro + probe allows setting gain independently for each channel
        % 3A or 3B data?
        % 3A metadata has field "typeEnabled" which was replaced
        % with "typeImEnabled" and "typeNiEnabled" in 3B.
        % The 3B imro table has an additional field for the
        % high pass filter enabled/disabled
        if isfield(meta,'typeEnabled')
          % 3A data
          C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
        else
          % 3B data
          C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d %*s', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
        end
        APgain = double(cell2mat(C(1)));
        LFgain = double(cell2mat(C(2)));
      else
        % get gain from  imChan0apGain, if present
        if isfield(meta,'imChan0apGain')
          APgain = APgain + str2num(meta.imChan0apGain);
          if acqCountList(2) > 0
            LFgain = LFgain + str2num(meta.imChan0lfGain);
          end
        elseif (probeType == 1110)
          % active UHD, for metadata lacking imChan0apGain, get gain from
          % imro table header
          currList = sscanf(meta.imroTbl, '(%d,%d,%d,%d,%d');
          APgain = APgain + currList(4);
          LFgain = LFgain + currList(5);
        elseif (probeType == 21) || (probeType == 24)
          % development NP 2.0; APGain = 80 for all AP
          % return 0 for LFgain (no LF channels)
          APgain = APgain + 80;
        elseif (probeType == 2013)
          % commercial NP 2.0; APGain = 80 for all AP
          APgain = APgain + 100;
        else
          fprintf('unknown gain, setting APgain to 1\n');
          APgain = APgain + 1;
        end
      end
    end % ChanGainsIM


    % =========================================================
    % Having acquired a block of raw nidq data using ReadBin(),
    % convert values to gain-corrected voltages. The conversion
    % is only applied to the saved-channel indices in chanList.
    % Remember saved-channel indices are in range [1:nSavedChans].
    % The dimensions of the dataArray remain unchanged. ChanList
    % examples:
    %
    %   [1:MN]      % all MN chans (MN from ChannelCountsNI)
    %   [2,6,20]    % just these three channels
    %
    function dataArray = GainCorrectNI(dataArray, chanList, meta)

      [MN,MA] = SGLX_readMeta.ChannelCountsNI(meta);
      fI2V = SGLX_readMeta.Int2Volts(meta);

      for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        conv = fI2V / SGLX_readMeta.ChanGainNI(j, MN, MA, meta);
        dataArray(j,:) = dataArray(j,:) * conv;
      end
    end % GainCorrectNI


    % =========================================================
    % Having acquired a block of raw imec data using ReadBin(),
    % convert values to gain-corrected voltages. The conversion
    % is only applied to the saved-channel indices in chanList.
    % Remember saved-channel indices are in range [1:nSavedChans].
    % The dimensions of the dataArray remain unchanged. ChanList
    % examples:
    %
    %   [1:AP]      % all AP chans (AP from ChannelCountsIM)
    %   [2,6,20]    % just these three channels
    %
    function dataArray = GainCorrectIM(dataArray, chanList, meta)

      % Look up gain with acquired channel ID
      chans = SGLX_readMeta.OriginalChans(meta);
      [APgain,LFgain] = SGLX_readMeta.ChanGainsIM(meta);
      nAP = length(APgain);
      nNu = nAP * 2;

      % Common conversion factor
      fI2V = SGLX_readMeta.Int2Volts(meta);

      for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        k = chans(j);       % acquisition index
        if k <= nAP
          conv = fI2V / APgain(k);
        elseif k <= nNu
          conv = fI2V / LFgain(k - nAP);
        else
          continue;
        end
        dataArray(j,:) = dataArray(j,:) * conv;
      end
    end % GainCorrectIM

    % =========================================================
    % Return array of survey bank times
    %
    %
    function bankTimes = svyBankTimes(meta)

      % Look up gain with acquired channel ID
      C = textscan(meta.svySBTT, '(%d %d %d %d', ...
        'EndOfLine', ')' );

      nBank = numel(C{1}) + 1;  % bank0/shank0 is at time = 0
      srate = SGLX_readMeta.SampRate(meta);

      bankTimes = zeros([nBank,4], "double");
      bankTimes(2:nBank,1) = double(C{1});
      bankTimes(2:nBank,2) = double(C{2});
      bankTimes(2:nBank,3) = double(C{3})/srate;
      bankTimes(2:nBank,4) = double(C{4})/srate;

    end % svyBankTimes

    % =========================================================
    % Write metadata file using values in meta structure
    %
    %
    function writeMeta(meta, newPath)

      % Write out metadata file. Tag order matches order of addition to
      % structure when read
      fmeta = fopen( newPath, 'w');
      tildeTags{1} = 'muxTbl';
      tildeTags{2} = 'imroTbl';
      tildeTags{3} = 'snsChanMap';
      tildeTags{4} = 'snsGeomMap';
      tildeTags{5} = 'snsShankMap';

      fn = fieldnames(meta);
      for i = 1:numel(fieldnames(meta))
        currTag = fn{i};
        tagFound = find(strcmp(tildeTags, currTag));
        if isempty(tagFound)
          currLine = sprintf('%s=%s',currTag,meta.(currTag));
          fprintf(fmeta,'%s\n',currLine);
        else
          currLine = sprintf('~%s=%s',currTag,meta.(currTag));
          fprintf(fmeta,'%s\n',currLine);
        end
      end

    end % writeMeta

    % ===========================================================
    % parse SGLX  imec filename with or without extension, return
    % runName,
    % gateStr, e.g. 'g0'
    % triggerStr, e.g. 't0' or 'tcat'
    % probeStr, e.g. 'imec0'
    % streamStr, e.g. 'ap' or 'lf'
    %
    %
    function [runName,gateStr,triggerStr,probeStr,streamStr] = parseFileName(fileName)

      % Remove extension, if present
      if endsWith(fileName, '.bin')
        fileName = fileName(1:length(fileName)-4);
      elseif endsWith(fileName, '.meta')
        fileName = fileName(1:length(fileName)-5);
      end

      % Find periods and underscores
      perPos = strfind(fileName,'.');
      usPos = strfind(fileName,'_');
      nPer = length(perPos);
      nUS = length(usPos);
      streamStr = fileName(perPos(nPer)+1:end);
      probeStr = fileName(perPos(nPer-1)+1:perPos(nPer)-1);
      triggerStr = fileName(usPos(nUS)+1:perPos(nPer-1)-1);
      gateStr = fileName(usPos(nUS-1)+1:usPos(nUS)-1);
      runName = fileName(1:usPos(nUS-1)-1);

    end % parseFileName

  end % SGLX_readMeta methods

end % SGLX_readMeta classdef
