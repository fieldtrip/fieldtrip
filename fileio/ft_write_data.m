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
%   'header'       = header structure that describes the data, see FT_READ_HEADER
%   'event'        = event structure that corresponds to the data, see FT_READ_EVENT
%   'chanindx'     = 1xN array, for selecting a subset of channels from header and data
%   'dataformat'   = string, see below
%   'append'       = boolean, not supported for all formats
%
% The supported dataformats for writing are
%   edf
%   gdf
%   anywave_ades
%   brainvision_eeg
%   neuralynx_ncs
%   neuralynx_sdma
%   plexon_nex
%   fcdc_matbin
%   fcdc_mysql
%   fcdc_buffer
%   flac, m4a, mp4, oga, ogg, wav (audio formats)
%   matlab
%   homer_nirs
%   snirf
%
% For EEG data, the input data is assumed to be scaled in microvolt.
% For NIRS data, the input data is assumed to represent optical densities.
%
% See also FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, FT_WRITE_EVENT

% Copyright (C) 2007-2021, Robert Oostenveld
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
nbits         = ft_getopt(varargin, 'nbits', 16); % for audio
chanindx      = ft_getopt(varargin, 'chanindx');
hdr           = ft_getopt(varargin, 'header');
evt           = ft_getopt(varargin, 'event');
dataformat    = ft_getopt(varargin, 'dataformat');

if isempty(dataformat)
  % only do the autodetection if the format was not specified
  dataformat = ft_filetype(filename);
end

if startsWith(dataformat, 'audio_')
  % support for  audio formats is implemented in a generic fashion
  dataformat = dataformat(7:end);
end

% convert 'yes' or 'no' string into boolean
append = istrue(append);

% ensure that the directory exists if we want to write to a file
if ~ismember(dataformat, {'empty', 'fcdc_global', 'fcdc_buffer', 'fcdc_mysql'})
  isdir_or_mkdir(fileparts(filename));
end

% determine the data size
[nchans, nsamples] = size(dat);

% ensure that the header is (reasonably) complete
if ~isfield(hdr, 'nChans')
  if isfield(hdr, 'label')
    hdr.nChans = length(hdr.label);
  else
    hdr.nChans = nchans;
  end
end

if ~isfield(hdr, 'label')
  hdr.label = arrayfun(@num2str, 1:hdr.nChans, 'UniformOutput', false)';
end

if ~isfield(hdr, 'chantype')
  % use a helper function which has some built in intelligence
  hdr.chantype = ft_chantype(hdr);
end

if ~isfield(hdr, 'chanunit')
  % use a helper function which has some built in intelligence
  hdr.chanunit = ft_chanunit(hdr);
end

if ~isempty(chanindx)
  % the header (and possibly the data) correspond to the original multichannel file
  hdr.label    = hdr.label(chanindx);
  hdr.chantype = hdr.chantype(chanindx);
  hdr.chanunit = hdr.chanunit(chanindx);
  hdr.nChans   = length(chanindx);
  if length(chanindx)==nchans
    % assume that the data already represents the desired subset of channels
  else
    % assume that the data corresponds to the original multichannel file
    dat = dat(chanindx,:);
  end
end

switch dataformat
  
  case 'empty'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % just pretend that we are writing the data, this is only for debugging
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [numC, numS] = size(dat);
    ft_info('Pretending to write %i samples from %i channels...\n',numS,numC);
    % Insert a small delay to make this more realitic for testing purposes
    % The time for writing to an actual location will differ and depend on
    % the amount of data
    pause(0.001);
    
  case 'fcdc_global'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store it in a global variable, this is only for debugging
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(evt)
      ft_error('writing events is not supported here, please see FT_WRITE_EVENT');
    end
    
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
    % write to a network transparent buffer for realtime analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(evt)
      ft_error('writing events is not supported here, please see FT_WRITE_EVENT');
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
        if contains(lasterr, 'Buffer size N must be an integer-valued scalar double.')
          % this happens if the MATLAB75/toolbox/signal/signal/buffer
          % function is used instead of the FieldTrip buffer
          ft_error('the FieldTrip buffer mex file was not found on your path, it should be in fieldtrip/fileio/private');
          
        elseif contains(lasterr, 'failed to create socket') && (strcmp(host, 'localhost') || strcmp(host, '127.0.0.1'))
          
          % start a local instance of the TCP server
          ft_warning('starting FieldTrip buffer on %s:%d', host, port);
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
      MAXNUMSAMPLE = 600000; % see buffer.h
      if size(dat,2)>MAXNUMSAMPLE
        ft_error('number of samples exceeds the size of the ring buffer');
      end
      % reformat the data into a buffer-compatible format
      packet.nchans    = size(dat,1);
      packet.nsamples  = size(dat,2);
      packet.data_type = find(strcmp(type, class(dat))) - 1; % zero-offset
      packet.bufsize   = numel(dat) * wordsize{strcmp(type, class(dat))};
      packet.buf       = dat;
      buffer('put_dat', packet, host, port);
    end
    
  case {'brainvision_eeg', 'brainvision_vhdr'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % combination of *.eeg and *.vhdr file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
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
      
      if isfield(old, 'event')
        event = old.event;
      else
        event = [];
      end
      % append the existing and the new events
      event = appendstruct(event, evt);
      
      save(headerfile, 'hdr', 'event', '-v6');
      
      % update the data file
      fid = fopen_or_error(datafile,'ab','ieee-le');
      fwrite(fid, dat, hdr.precision);
      fclose(fid);
      
    else
      hdr.nSamples = nsamples;
      hdr.nTrials  = 1;
      if ~isfield(hdr, 'precision')
        hdr.precision = 'double';
      end
      
      % rename the variable name for the new events
      event = evt;
      
      % write the header and events to the file
      save(headerfile, 'hdr', 'event', '-v6');
      
      % write the data file
      fid = fopen_or_error(datafile,'wb','ieee-le');
      fwrite(fid, dat, hdr.precision);
      fclose(fid);
    end
    
    
  case 'fcdc_mysql'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write to a MySQL server listening somewhere else on the network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(evt)
      ft_error('writing events is not supported here, please see FT_WRITE_EVENT');
    end
    
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
          ft_warning(lasterr);
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
            ft_warning(lasterr);
          end
          % insert the structure into the database
          db_insert('fieldtrip.data', s);
        end
      end
      
    else
      ft_error('you should specify either the header or the data when writing to a MySQL database');
    end
    
  case 'matlab'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plain MATLAB file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file '.mat']);
    
    if  append &&  exist(filename, 'file')
      % read the previous header and data from MATLAB file
      prev = load(filename);
      % do a sanity chjeck to ensure that the file content is consistent with the new data
      if ~isempty(hdr) && ~isequal(hdr, prev.hdr)
        ft_error('inconsistent header');
      end
      
    elseif append && ~exist(filename, 'file')
      % file does not yet exist, which is not a problem
      prev = [];
      
    elseif ~append &&  exist(filename, 'file')
      % file already exists, delete it and make a new one further down
      ft_warning('deleting existing file ''%s''', filename);
      delete(filename);
      prev = [];
      
    elseif ~append && ~exist(filename, 'file')
      % file does not yet exist, which is not a problem
      prev = [];
    end
    
    if isfield(prev, 'dat')
      % append the new data to the previous data from from the MATLAB file
      dat = cat(2, prev.dat, dat);
    end
    if isfield(prev, 'event')
      % append the new events to the previous events from from the MATLAB file
      event = cat(2, prev.event, evt);
    else
      % rename the variable name for the new events
      event = evt;
    end
    
    % write the data, header and events to the file
    save(filename, 'dat', 'hdr', 'event');
    
  case 'mff'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MFF files using Phillips plugin
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ft_hastoolbox('mffmatlabio', 1);
    mff_fileio_write(filename, hdr, dat, evt);
    
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
    if ~isempty(evt)
      ft_error('writing events is not supported');
    end
    
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
    
    if ~isfolder(dirname)
      mkdir(dirname);
    end
    
    % open and write to the output files, one for each selected channel
    fid = zeros(hdr.nChans,1);
    for j=1:hdr.nChans
      
      if append==false
        fid(j) = fopen_or_error(filename{j}, 'wb', 'ieee-le'); % open the file
        magic = format{j};                               % this used to be the channel name
        magic((end+1):8) = 0;                            % pad with zeros
        magic(8) = downscale(j);                         % number of bits to shift
        fwrite(fid(j), magic(1:8));                      % write the 8-byte file header
      else
        fid(j) = fopen_or_error(filename{j}, 'ab', 'ieee-le');    % open the file for appending
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
            ft_error('unsupported format conversion');
        end
      end
      
      % apply the scaling, this corresponds to bit shifting
      buf = buf ./ (2^downscale(j));
      
      % write the segment of data to the output file
      fwrite(fid(j), buf, format{j}, 'ieee-le');
      
      fclose(fid(j));
    end % for each channel
    
  case {'flac' 'm4a' 'mp4' 'oga' 'ogg' 'wav'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This writes data Y to a Windows WAVE file specified by the file name
    % WAVEFILE, with a sample rate of FS Hz and with NBITS number of bits.
    % NBITS must be 8, 16, 24, or 32.  For NBITS < 32, amplitude values
    % outside the range [-1,+1] are clipped
    %
    % Supported extensions for AUDIOWRITE are:
    %   .flac
    % 	.m4a
    % 	.mp4
    % 	.oga
    % 	.ogg
    % 	.wav
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not supported for this data format');
    end
    if ~isempty(evt)
      ft_error('writing events is not supported');
    end
    
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file '.' dataformat]);
    
    switch dataformat
      case {'m4a' 'mp4' 'oga' 'ogg'}
        options = {};
      otherwise
        options = {'BitsPerSample', nbits};
    end % switch
    
    audiowrite(filename, dat', hdr.Fs, options{:});
    
  case 'plexon_nex'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single or mulitple channel Plexon NEX file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end
    if ~isempty(evt)
      ft_error('writing events is not supported');
    end
    if nchans~=1
      ft_error('only supported for single-channel data');
    end
    
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.nex']);
    
    % construct a NEX structure with  the required parts of the header
    nex.hdr.VarHeader.Type       = 5; % continuous
    nex.hdr.VarHeader.Name       = hdr.label{1};
    nex.hdr.VarHeader.WFrequency = hdr.Fs;
    if isfield(hdr, 'FirstTimeStamp')
      nex.hdr.FileHeader.Frequency = hdr.Fs * hdr.TimeStampPerSample;
      nex.var.ts = hdr.FirstTimeStamp;
    else
      ft_warning('no timestamp information available');
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
      ft_error('appending data is not yet supported for this data format');
    end
    if ~isempty(evt)
      ft_error('writing events is not supported');
    end
    if nchans>1
      ft_error('only supported for single-channel data');
    end
    
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.ncs']);
    
    LABEL      = hdr.label{1};  % single channel
    ADCHANNEL  = -1;            % unknown
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
    
    if false
      % the following code snippet can be used for testing
      ncs2 = read_neuralynx_ncs(filename, 1, inf);
    end
    
  case 'gdf'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiple channel GDF file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end
    if ~isempty(evt)
      ft_error('writing events is not supported');
    end
    
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.gdf']);
    
    write_gdf(filename, hdr, dat);
    
  case 'edf'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiple channel European Data Format file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end
    if ~isempty(evt)
      ft_error('writing events is not supported');
    end
    
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.edf']);
    
    write_edf(filename, hdr, dat);
    
  case 'anywave_ades'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see http://meg.univ-amu.fr/wiki/AnyWave:ADES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end
    if ~isempty(evt)
      ft_error('writing events is not supported');
    end
    
    dattype = unique(hdr.chantype);
    datunit = cell(size(dattype));
    for i=1:numel(dattype)
      unit = hdr.chanunit(strcmp(hdr.chantype, dattype{i}));
      if ~all(strcmp(unit, unit{1}))
        ft_error('channels of the same type with different units are not supported');
      end
      datunit{i} = unit{1};
    end
    
    % only change these after checking channel types and units
    chantype = adestype(hdr.chantype);
    dattype  = adestype(dattype);
    
    % ensure that all channels have the right scaling
    for i=1:size(dat,1)
      switch chantype{i}
        case 'MEG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'T');
        case 'Reference'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'T');
        case 'GRAD'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'T/m');
        case 'EEG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'V');
        case 'SEEG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'V');
        case 'EMG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'V');
        case 'ECG'
          dat(i,:) = dat(i,:) * ft_scalingfactor(hdr.chanunit{i}, 'V');
        otherwise
          % FIXME I am not sure what scaling to apply
      end
    end
    
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, f); % without extension
    mat2ades(dat, filename, hdr.Fs, hdr.label, chantype, dattype, datunit);
    
  case 'homer_nirs'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % https://www.nitrc.org/plugins/mwiki/index.php/homer2:Homer_Input_Files#NIRS_data_file_format
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end
    
    % convert the input arguments into a FieldTrip raw data structure
    data = [];
    data.hdr = hdr;
    data.trial{1} = dat;
    data.time{1} = ((1:hdr.nSamples*hdr.nTrials)-1)/hdr.Fs;
    data.label = hdr.label;
    data.sampleinfo = [1 size(dat,2)];
    
    % convert the raw data structure to Homer format
    nirs = fieldtrip2homer(data, 'event', evt);
    
    % Homer files are MATLAB files in disguise
    % see https://www.nitrc.org/plugins/mwiki/index.php/homer2:Homer_Input_Files#NIRS_data_file_format
    save(filename, '-struct', 'nirs'); % save the fields as individual variables in the file
    
  case 'snirf'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % https://github.com/fNIRS/snirf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if append
      ft_error('appending data is not yet supported for this data format');
    end
    
    % this uses the SNIRF reading functions from the Homer3 toolbox
    ft_hastoolbox('homer3', 1);
    
    % construct a time axis that matches the data, it starts at 0 seconds
    time = ((1:hdr.nSamples*hdr.nTrials)-1)/hdr.Fs;
    
    % divide data in nirs channels, stimulus channels and auxillary channels
    seldat  = startsWith(hdr.chantype, 'nirs');
    selstim = strcmp(hdr.chantype, 'stimulus');
    selaux  = ~seldat & ~selstim;
    
    % create empty SnirfClass
    snirf = SnirfClass();
    
    % collect information for creation of snirf file
    source_idx   = find(contains(hdr.opto.optotype, {'transmitter', 'source'}));
    detector_idx = find(contains(hdr.opto.optotype, {'receiver', 'detector'}));
    tra = hdr.opto.tra;
    tra_t = tra'; % transpose tra matrix to get indices of wavelength by row (thus by channel)
    wl_idx = find(tra_t>0);
    all_wavelengths = hdr.opto.wavelength(tra_t(wl_idx));
    split = nanmedian(all_wavelengths);
    WL1.values = all_wavelengths(all_wavelengths<split);
    WL2.values = all_wavelengths(all_wavelengths>split);
    WL1.nominal = round(median(WL1.values),-1);
    WL2.nominal = round(median(WL2.values),-1);
    ft_warning('assuming that the nominal wavelengths are %d and %d nm', WL1.nominal, WL2.nominal)
    
    % metaDataTags
    snirf.metaDataTags(1).tags.LengthUnit = hdr.opto.unit;
    snirf.metaDataTags(1).tags.TimeUnit = 's';
    snirf.metaDataTags(1).tags.FrequencyUnit = 'Hz';
    
    % data
    snirf.data(1).dataTimeSeries = dat(seldat,:)'; % <number of time points> x <number of channels>
    snirf.data(1).time = time';                    % <number of time points x 1> (can also be  represented as <start time x sample time spacing>
    
    % measurementList
    for i=1:size(tra,1)
      source = find(tra(i,:)>0);
      detector = find(tra(i,:)<0);
      snirf.data.measurementList(i).sourceIndex   = find(source_idx==source);
      snirf.data.measurementList(i).detectorIndex = find(detector_idx==detector);
      %     snirf.data.measurementList(i).wavelengthActual = all_wavelengths(i); % this is not yet supported by the snirf toolbox
      if any(all_wavelengths(i)==WL1.values)
        snirf.data.measurementList(i).wavelengthIndex = 1;
      else
        snirf.data.measurementList(i).wavelengthIndex = 2;
      end
      snirf.data.measurementList(i).dataType = 99999;
      snirf.data.measurementList(i).dataTypeLabel = 'dOD';
    end
    ft_warning('assuming that the input data represents (change in) optical densities')
    
    % sort the channels according to wavelengths, because this is the way that homer handles data
    [dum, idx] = sort(([snirf.data.measurementList(:).wavelengthIndex]));
    % update the data accordingly
    snirf.data.measurementList = snirf.data.measurementList(idx);
    snirf.data.dataTimeSeries=snirf.data.dataTimeSeries(:, idx);
    
    % stim
    if ~isempty(evt)
      % select events with a string value, the type will be a string
      sel = cellfun(@ischar, {evt.value});
      if any(sel)
        evt_string = evt(sel);
        evt_names = unique({evt_string(:).value}); % if the values are strings, this propably contains the event names
        for i=1:length(evt_names)
          snirf.stim(i).name = evt_names{i};
          evt_idx = find(strcmp({evt_string(:).value}, evt_names{i}));
          starttime = ([evt_string(evt_idx).sample]-1)/hdr.Fs;
          duration = [evt_string(evt_idx).duration]/hdr.Fs;
          if isempty(duration)
            duration = zeros(1, length(starttime));
          end
          value = ones(1,length(evt_idx));
          snirf.stim(i).data = [starttime' duration' value'];
        end
      end
      % select events with a numeric value, the type will be a string
      sel = cellfun(@isnumeric, {evt.value});
      if any(sel)
        evt_numeric = evt(sel);
        evt_names = unique({evt_numeric(:).type});
        for i=1:length(evt_names)
          snirf.stim(i).name = evt_names{i};
          evt_idx = find(strcmp({evt_numeric(:).type}, evt_names{i}));
          starttime = ([evt_numeric(evt_idx).sample]-1)/hdr.Fs;
          duration = [evt_numeric(evt_idx).duration]/hdr.Fs;
          if isempty(duration)
            duration = zeros(1, length(starttime));
          end
          value = [evt_numeric(evt_idx).value];
          if isempty(value)
            value = ones(1, length(starttime));
          end
          snirf.stim(i).data = [starttime' duration' value'];
        end
      end
    end
    
    % probe
    snirf.probe(1).wavelengths = [WL1.nominal WL2.nominal];
    if all(hdr.opto.optopos(:,3)==0)
      snirf.probe(1).sourcePos2D    = hdr.opto.optopos(source_idx, 1:2);
      snirf.probe(1).detectorPos2D  = hdr.opto.optopos(detector_idx, 1:2);
    else
      snirf.probe(1).sourcePos3D    = hdr.opto.optopos(source_idx, 1:3);
      snirf.probe(1).detectorPos3D  = hdr.opto.optopos(detector_idx, 1:3);
      layoutpos = getorthoviewpos(hdr.opto.optopos, 'ras', 'superior');
      snirf.probe(1).sourcePos2D    = layoutpos(source_idx, 1:2);
      snirf.probe(1).detectorPos2D  = layoutpos(detector_idx, 1:2);
    end
    snirf.probe(1).sourceLabels     = hdr.opto.optolabel(source_idx);
    snirf.probe(1).detectorLabels   = hdr.opto.optolabel(detector_idx);
    
    % aux
    if sum(selaux)~=0
      auxdata = dat(selaux,:);
      auxlabels = hdr.label(selaux);
      for i=1:sum(selaux)
        snirf.aux(i).name = auxlabels{i}; % check if correct format
        snirf.aux(i).dataTimeSeries = auxdata(i,:)';
        snirf.aux(i).time = time';
      end
    end
    
    % save .snirf file
    snirf.Save(filename)
    
  otherwise
    ft_error('unsupported data format');
end % switch dataformat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function type = adestype(type)
for i=1:numel(type)
  switch lower(type{i})
    case 'meggrad'
      type{i} = 'MEG'; % this is for CTF and BTi/4D
    case 'megmag'
      type{i} = 'MEG'; % this is for Neuromag and BTi/4D
    case 'megplanar'
      type{i} = 'GRAD'; % this is for Neuromag
    case {'refmag' 'refgrad'}
      type{i} = 'Reference';
    case 'eeg'
      type{i} = 'EEG';
    case {'seeg' 'ecog' 'ieeg'}
      type{i} = 'SEEG'; % all intracranial channels
    case 'ecg'
      type{i} = 'ECG';
    case 'emg'
      type{i} = 'EMG';
    case 'trigger'
      type{i} = 'Trigger';
    case 'source'
      type{i} = 'Source'; % virtual channel
    otherwise
      type{i} = 'Other';
  end
end
