function ft_write_data(filename, dat, varargin)

% FT_WRITE_DATA exports electrophysiological data such as EEG to a file. 
%
% Use as
%   ft_write_data(filename, dat, ...)
%
% The specified filename can contain the filename extension. If it has no filename
% extension not, it will be added automatically.
%
% Additional options should be specified in key-value pairs and can be
%   'header'         header structure that describes the data, see FT_READ_HEADER
%   'dataformat'     string, see below
%   'append'         boolean, not supported for all formats
%   'chanindx'       1xN array
%
% The supported dataformats are
%   edf
%   gdf
%   brainvision_eeg
%   neuralynx_ncs
%   neuralynx_sdma
%   plexon_nex
%   riff_wave
%   fcdc_matbin
%   fcdc_mysql
%   fcdc_buffer
%   matlab
%
% For EEG data formats, the input data is assumed to be scaled in microvolt.
%
% See also FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, FT_WRITE_EVENT

% Copyright (C) 2007-2014, Robert Oostenveld
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

global data_queue    % for fcdc_global
global header_queue  % for fcdc_global
global db_blob       % for fcdc_mysql
if isempty(db_blob)
  db_blob = 0;
end

% get the options
append        = ft_getopt(varargin, 'append', false);
nbits         = ft_getopt(varargin, 'nbits', 16); % for riff_wave
chanindx      = ft_getopt(varargin, 'chanindx');
hdr           = ft_getopt(varargin, 'header');
evt           = ft_getopt(varargin, 'event');
dataformat    = ft_getopt(varargin, 'dataformat');

if isempty(dataformat)
  % only do the autodetection if the format was not specified
  dataformat = ft_filetype(filename);
end

% convert 'yes' or 'no' string into boolean
append = istrue(append);

% determine the data size
[nchans, nsamples] = size(dat);

switch dataformat
  
  case 'empty'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % just pretend that we are writing the data, this is only for debugging
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [numC, numS] = size(dat);
    fprintf(1,'Pretending to write %i samples from %i channels...\n',numS,numC);
    % Insert a small delay to make this more realitic for testing purposes
    % The time for writing to an actual location will differ and depend on
    % the amount of data
    pause(0.001);
    
  case 'fcdc_global'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store it in a global variable, this is only for debugging
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(hdr)
      header_queue = hdr;
    end
    if isempty(data_queue) || ~append
      data_queue = dat;
    else
      data_queue = cat(2, data_queue, dat);
    end
    
  case 'fcdc_buffer'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % network transparent buffer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(chanindx)
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    
    [host, port] = filetype_check_uri(filename);
    
    type = {
      'char'
      'uint8'
      'uint16'
      'uint32'
      'uint64'
      'int8'
      'int16'
      'int32'
      'int64'
      'single'
      'double'
      };
    
    wordsize = {
      1 % 'char'
      1 % 'uint8'
      2 % 'uint16'
      4 % 'uint32'
      8 % 'uint64'
      1 % 'int8'
      2 % 'int16'
      4 % 'int32'
      8 % 'int64'
      4 % 'single'
      8 % 'double'
      };
    
    % this should only be done the first time
    if ~append && ~isempty(hdr)
      % reformat the header into a buffer-compatible format
      packet.fsample   = hdr.Fs;
      packet.nchans    = hdr.nChans;
      packet.nsamples  = 0;
      packet.nevents   = 0;
      packet.data_type = find(strcmp(type, class(dat))) - 1; % zero-offset
      if isfield(hdr,'label') && iscell(hdr.label)
        packet.channel_names = hdr.label;
      end
      if isfield(hdr,'siemensap')
        if isa(hdr.siemensap, 'uint8')
          packet.siemensap = hdr.siemensap;
        else
          %          try
          %            packet.siemensap = matlab2sap(hdr.siemensap);
          %          catch
          warning 'Ignoring field "siemensap"';
          %          end
        end
      end
      if isfield(hdr,'nifti_1')
        if isa(hdr.nifti_1, 'uint8')
          packet.nifti_1 = hdr.nifti_1;
        else
          try
            packet.nifti_1 = encode_nifti1(hdr.nifti_1);
          catch
            warning 'Ignoring field "nifti_1"';
          end
        end
      end
      if isfield(hdr,'ctf_res4')
        if isa(hdr.ctf_res4, 'uint8')
          packet.ctf_res4 = hdr.ctf_res4;
        else
          warning 'Ignoring non-uint8 field "ctf_res4"';
        end
      end
      
      % try to put_hdr and initialize if necessary
      try
        % try writing the packet
        buffer('put_hdr', packet, host, port);
        
      catch
        if ~isempty(strfind(lasterr, 'Buffer size N must be an integer-valued scalar double.'))
          % this happens if the MATLAB75/toolbox/signal/signal/buffer
          % function is used instead of the FieldTrip buffer
          error('the FieldTrip buffer mex file was not found on your path, it should be in fieldtrip/fileio/private');
          
        elseif ~isempty(strfind(lasterr, 'failed to create socket')) && (strcmp(host, 'localhost') || strcmp(host, '127.0.0.1'))
          
          % start a local instance of the TCP server
          warning('starting FieldTrip buffer on %s:%d', host, port);
          buffer('tcpserver', 'init', host, port);
          pause(1);
          
          % rewrite the packet until success
          success = false;
          while ~success
            try
              % it may take some time before the TCP server is fully initialized
              % try writing the packet again
              buffer('put_hdr', packet, host, port);
              success = true;
            catch
              success = false;
            end
          end
        end % if strfind...
        
      end % try
      
    end % writing header
    
    if ~isempty(dat)
      max_nsamples = 32556;
      if size(dat,2)>max_nsamples
        % FIXME this is a hack to split large writes into multiple smaller writes
        % this is to work around a problem observed in the neuralynx proxy
        % when sampling 32 channels at 32KHz
        begsample = 1;
        while begsample<=size(dat,2)
          endsample = begsample - 1 + max_nsamples;
          endsample = min(endsample, size(dat,2));
          % if append is already one of the arguments, remove it from varargin
          indx = find(strcmp(varargin, 'append')); % find the "append" key
          if ~isempty(indx)
            indx = [indx indx+1];                  % remove the key and the value
            varargin(indx) = [];
          end
          ft_write_data(filename, dat(:,begsample:endsample), varargin{:}, 'append', false);
          begsample = endsample + 1;
        end
      else
        % FIXME this is the normal code, which will also be used recursively
        % reformat the data into a buffer-compatible format
        packet.nchans    = size(dat,1);
        packet.nsamples  = size(dat,2);
        packet.data_type = find(strcmp(type, class(dat))) - 1; % zero-offset
        packet.bufsize   = numel(dat) * wordsize{find(strcmp(type, class(dat)))};
        packet.buf       = dat;
        buffer('put_dat', packet, host, port);
      end % if data larger than chuncksize
    end
    
  case 'brainvision_eeg'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % combination of *.eeg and *.vhdr file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      error('appending data is not yet supported for this data format');
    end
    
    if nchans~=hdr.nChans && length(chanindx)==nchans
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    % the header should at least contain the following fields
    %   hdr.label
    %   hdr.nChans
    %   hdr.Fs
    write_brainvision_eeg(filename, hdr, dat, evt);
    
  case 'fcdc_matbin'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiplexed data in a *.bin file (ieee-le, 64 bit floating point values),
    % accompanied by a MATLAB V6 file containing the header
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.mat']);
    datafile   = fullfile(path, [file '.bin']);
    
    if append && exist(headerfile, 'file') && exist(datafile, 'file')
      % read the existing header and perform a sanity check
      old = load(headerfile);
      assert(old.hdr.nChans==size(dat,1));
      
      % update the existing header
      hdr          = old.hdr;
      hdr.nSamples = hdr.nSamples + nsamples;
      
      % there are no new events
      if isfield(old, 'event')
        event = old.event;
      else
        event = [];
      end
      
      save(headerfile, 'hdr', 'event', '-v6');

      % update the data file
      [fid,message] = fopen(datafile,'ab','ieee-le');
      fwrite(fid, dat, hdr.precision);
      fclose(fid);
      
    else
      hdr.nSamples = nsamples;
      hdr.nTrials  = 1;
      if nchans~=hdr.nChans && length(chanindx)==nchans
        % assume that the header corresponds to the original multichannel
        % file and that the data represents a subset of channels
        hdr.label     = hdr.label(chanindx);
        hdr.nChans    = length(chanindx);
      end
      if ~isfield(hdr, 'precision')
        hdr.precision = 'double';
      end
      % there are no events
      event = [];
      % write the header file
      save(headerfile, 'hdr', 'event', '-v6');
      
      % write the data file
      [fid,message] = fopen(datafile,'wb','ieee-le');
      fwrite(fid, dat, hdr.precision);
      fclose(fid);
    end
    
    
  case 'fcdc_mysql'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write to a MySQL server listening somewhere else on the network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check that the required low-level toolbox is available
    ft_hastoolbox('mysql', 1);
    db_open(filename);
    if ~isempty(hdr) && isempty(dat)
      % insert the header information into the database
      if db_blob
        % insert the structure into the database table as a binary blob
        db_insert_blob('fieldtrip.header', 'msg', hdr);
      else
        % make a structure with the same elements as the fields in the database table
        s             = struct;
        s.Fs          = hdr.Fs;           % sampling frequency
        s.nChans      = hdr.nChans;       % number of channels
        s.nSamples    = hdr.nSamples;     % number of samples per trial
        s.nSamplesPre = hdr.nSamplesPre;  % number of pre-trigger samples in each trial
        s.nTrials     = hdr.nTrials;      % number of trials
        s.label       = mxSerialize(hdr.label);
        try
          s.msg = mxSerialize(hdr);
        catch
          warning(lasterr);
        end
        db_insert('fieldtrip.header', s);
      end
      
    elseif isempty(hdr) && ~isempty(dat)
      dim = size(dat);
      if numel(dim)==2
        % ensure that the data dimensions correspond to ntrials X nchans X samples
        dim = [1 dim];
        dat = reshape(dat, dim);
      end
      ntrials = dim(1);
      for i=1:ntrials
        if db_blob
          % insert the data into the database table as a binary blob
          db_insert_blob('fieldtrip.data', 'msg', reshape(dat(i,:,:), dim(2:end)));
        else
          % create a structure with the same fields as the database table
          s = struct;
          s.nChans   = dim(2);
          s.nSamples = dim(3);
          try
            s.data = mxSerialize(reshape(dat(i,:,:), dim(2:end)));
          catch
            warning(lasterr);
          end
          % insert the structure into the database
          db_insert('fieldtrip.data', s);
        end
      end
      
    else
      error('you should specify either the header or the data when writing to a MySQL database');
    end
    
  case 'matlab'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plain MATLAB file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file '.mat']);
    if      append &&  exist(filename, 'file')
      % read the previous header and data from MATLAB file
      prev = load(filename);
      if ~isempty(hdr) && ~isequal(hdr, prev.hdr)
        error('inconsistent header');
      else
        % append the new data to that from the MATLAB file
        dat = cat(2, prev.dat, dat);
      end
    elseif  append && ~exist(filename, 'file')
      % file does not yet exist, which is not a problem
    elseif ~append &&  exist(filename, 'file')
      warning('deleting existing file ''%s''', filename);
      delete(filename);
    elseif ~append && ~exist(filename, 'file')
      % file does not yet exist, which is not a problem
    end
    save(filename, 'dat', 'hdr');
    
  case 'neuralynx_sdma'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The first version of this file format contained in the first 8 bytes the
    % channel label as string. Subsequently it contained 32 bit integer values.
    %
    % The second version of this file format starts with 8 bytes describing (as
    % a space-padded string) the data type. The channel label is contained in
    % the filename as dataset.chanlabel.bin.
    %
    % The third version of this file format starts with 7 bytes describing (as
    % a zero-padded string) the data type, followed by the 8th byte which
    % describes the downscaling for the 8 and 16 bit integer representations.
    % The downscaling itself is represented as uint8 and should be interpreted as
    % the number of bits to shift. The channel label is contained in the
    % filename as dataset.chanlabel.bin.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    statuschannel = {
      'stx'
      'pid'
      'siz'
      'tsh'
      'tsl'
      'cpu'
      'ttl'
      'x01'
      'x02'
      'x03'
      'x04'
      'x05'
      'x06'
      'x07'
      'x08'
      'x09'
      'x10'
      'crc'
      };
    
    dirname = filename;
    clear filename
    [path, file] = fileparts(dirname);
    for i=1:hdr.nChans
      downscale(i) = 0;
      if ~isempty(strmatch(hdr.label{i}, statuschannel))
        format{i} = 'uint32';
      else
        format{i} = 'int32';
      end
      filename{i} = fullfile(dirname, [file '.' hdr.label{i} '.bin']);
    end
    
    if ~isdir(dirname)
      mkdir(dirname);
    end
    
    % open and write to the output files, one for each selected channel
    fid = zeros(hdr.nChans,1);
    for j=1:hdr.nChans
      
      if append==false
        fid(j) = fopen(filename{j}, 'wb', 'ieee-le');   % open the file
        magic = format{j};                              % this used to be the channel name
        magic((end+1):8) = 0;                           % pad with zeros
        magic(8) = downscale(j);                        % number of bits to shift
        fwrite(fid(j), magic(1:8));                     % write the 8-byte file header
      else
        fid(j) = fopen(filename{j}, 'ab', 'ieee-le');   % open the file for appending
      end % if append
      
      % convert the data into the correct class
      buf = dat(j,:);
      if ~strcmp(class(buf), format{j})
        switch format{j}
          case 'int16'
            buf = int16(buf);
          case 'int32'
            buf = int32(buf);
          case 'single'
            buf = single(buf);
          case 'double'
            buf = double(buf);
          case 'uint32'
            buf = uint32(buf);
          otherwise
            error('unsupported format conversion');
        end
      end
      
      % apply the scaling, this corresponds to bit shifting
      buf = buf ./ (2^downscale(j));
      
      % write the segment of data to the output file
      fwrite(fid(j), buf, format{j}, 'ieee-le');
      
      fclose(fid(j));
    end % for each channel
    
  case 'riff_wave'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     This writes data Y to a Windows WAVE file specified by the file name
    %     WAVEFILE, with a sample rate of FS Hz and with NBITS number of bits.
    %     NBITS must be 8, 16, 24, or 32.  For NBITS < 32, amplitude values
    %     outside the range [-1,+1] are clipped
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      error('appending data is not yet supported for this data format');
    end
    
    if nchans~=hdr.nChans && length(chanindx)==nchans
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    if nchans~=1
      error('this format only supports single channel continuous data');
    end
    audiowrite(filename, dat, hdr.Fs, 'BitsPerSample', nbits);
    
  case 'plexon_nex'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single or mulitple channel Plexon NEX file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      error('appending data is not yet supported for this data format');
    end
    
    [path, file] = fileparts(filename);
    filename = fullfile(path, [file, '.nex']);
    if nchans~=1
      error('only supported for single-channel data');
    end
    % construct a NEX structure with  the required parts of the header
    nex.hdr.VarHeader.Type       = 5; % continuous
    nex.hdr.VarHeader.Name       = hdr.label{1};
    nex.hdr.VarHeader.WFrequency = hdr.Fs;
    if isfield(hdr, 'FirstTimeStamp')
      nex.hdr.FileHeader.Frequency = hdr.Fs * hdr.TimeStampPerSample;
      nex.var.ts = hdr.FirstTimeStamp;
    else
      warning('no timestamp information available');
      nex.hdr.FileHeader.Frequency  = nan;
      nex.var.ts = nan;
    end
    nex.var.indx = 0;
    nex.var.dat  = dat;
    
    write_plexon_nex(filename, nex);
    
    if 0
      % the following code snippet can be used for testing
      [nex2.var, nex2.hdr] = read_plexon_nex(filename, 'channel', 1);
    end
    
  case 'neuralynx_ncs'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single channel Neuralynx NCS file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      error('appending data is not yet supported for this data format');
    end
    
    if nchans>1
      error('only supported for single-channel data');
    end
    
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.ncs']);
    
    if nchans~=hdr.nChans && length(chanindx)==nchans
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      % WARNING the AD channel index assumes that the data was read from a DMA or SDMA file
      % the first 17 channels contain status info, this number is zero-offset
      ADCHANNEL  = chanindx - 17 - 1;
      LABEL      = hdr.label{chanindx};
    elseif hdr.nChans==1
      ADCHANNEL  = -1;            % unknown
      LABEL      = hdr.label{1};  % single channel
    else
      error('cannot determine channel label');
    end
    
    FSAMPLE    = hdr.Fs;
    RECORDNSMP = 512;
    RECORDSIZE = 1044;
    
    % cut the downsampled LFP data into record-size pieces
    nrecords = ceil(nsamples/RECORDNSMP);
    fprintf('construct ncs with %d records\n', nrecords);
    
    % construct a ncs structure with all header details and the data in it
    ncs                = [];
    ncs.NumValidSamp   = ones(1,nrecords) * RECORDNSMP;   % except for the last block
    ncs.ChanNumber     = ones(1,nrecords) * ADCHANNEL;
    ncs.SampFreq       = ones(1,nrecords) * FSAMPLE;
    ncs.TimeStamp      = zeros(1,nrecords,'uint64');
    
    if rem(nsamples, RECORDNSMP)>0
      % the data length is not an integer number of records, pad the last record with zeros
      dat = cat(2, dat, zeros(nchans, nrecords*RECORDNSMP-nsamples));
      ncs.NumValidSamp(end) = rem(nsamples, RECORDNSMP);
    end
    
    ncs.dat = reshape(dat, RECORDNSMP, nrecords);
    
    for i=1:nrecords
      % timestamps should be 64 bit unsigned integers
      ncs.TimeStamp(i) = uint64(hdr.FirstTimeStamp) + uint64((i-1)*RECORDNSMP*hdr.TimeStampPerSample);
    end
    
    % add the elements that will go into the ascii header
    ncs.hdr.CheetahRev            = '4.23.0';
    ncs.hdr.NLX_Base_Class_Type   = 'CscAcqEnt';
    ncs.hdr.NLX_Base_Class_Name   = LABEL;
    ncs.hdr.RecordSize            = RECORDSIZE;
    ncs.hdr.ADChannel             = ADCHANNEL;
    ncs.hdr.SamplingFrequency     = FSAMPLE;
    
    % write it to a file
    fprintf('writing to %s\n', filename);
    write_neuralynx_ncs(filename, ncs);
    
    if 0
      % the following code snippet can be used for testing
      ncs2 = read_neuralynx_ncs(filename, 1, inf);
    end
    
  case 'gdf'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiple channel GDF file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      error('appending data is not yet supported for this data format');
    end
    if ~isempty(chanindx)
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    write_gdf(filename, hdr, dat);
    
  case 'edf'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiple channel European Data Format file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      error('appending data is not yet supported for this data format');
    end
    if ~isempty(chanindx)
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    write_edf(filename, hdr, dat);
    
  otherwise
    error('unsupported data format');
end % switch dataformat
