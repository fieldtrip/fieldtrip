function [hdr] = ft_read_header(filename, varargin)

% FT_READ_HEADER reads header information from a variety of EEG, MEG and other time
% series data files and represents the header information in a common
% data-independent structure. The supported formats are listed below.
%
% Use as
%   hdr = ft_read_header(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'headerformat'   = string
%   'fallback'       = can be empty or 'biosig' (default = [])
%   'checkmaxfilter' = boolean, whether to check that maxfilter has been correctly applied (default = true)
%   'chanindx'       = list with channel indices in case of different sampling frequencies (only for EDF)
%   'chantype'       = string or cell-array with strings, channel types to be read (only for NeuroOmega and BlackRock)
%   'coordsys'       = string, 'head' or 'dewar' (default = 'head')
%   'headerformat'   = name of a MATLAB function that takes the filename as input (default is automatic)
%   'password'       = password structure for encrypted data set (only for mayo_mef30 and mayo_mef21)
%   'readbids'       = string, 'yes', no', or 'ifmakessense', whether to read information from the BIDS sidecar files (default = 'ifmakessense')
%
% This returns a header structure with the following fields
%   hdr.Fs          = sampling frequency
%   hdr.nChans      = number of channels
%   hdr.nSamples    = number of samples per trial
%   hdr.nSamplesPre = number of pre-trigger samples in each trial
%   hdr.nTrials     = number of trials
%   hdr.label       = Nx1 cell-array with the label of each channel
%   hdr.chantype    = Nx1 cell-array with the channel type, see FT_CHANTYPE
%   hdr.chanunit    = Nx1 cell-array with the physical units, see FT_CHANUNIT
%
% For continuously recorded data, this will return nSamplesPre=0 and nTrials=1.
%
% For some data formats that are recorded on animal electrophysiology
% systems (e.g. Neuralynx, Plexon), the following optional fields are
% returned, which allows for relating the timing of spike and LFP data
%   hdr.FirstTimeStamp      number, represented as 32-bit or 64-bit unsigned integer
%   hdr.TimeStampPerSample  number, represented in double precision
%
% Depending on the file format, additional header information can be
% returned in the hdr.orig subfield.
%
% To use an external reading function, you can specify an external function as the
% 'headerformat' option. This function should take the filename as input argument.
% Please check the code of this function for details, and search for BIDS_TSV as
% example.
%
% The following MEG dataformats are supported
%   CTF (*.ds, *.res4, *.meg4)
%   Neuromag/Elekta/Megin (*.fif)
%   BTi/4D (*.m4d, *.pdf, *.xyz)
%   Yokogawa/Ricoh (*.ave, *.con, *.raw)
%   NetMEG (*.nc)
%   ITAB - Chieti (*.mhd)
%   Tristan Babysquid (*.fif)
%
% The following EEG dataformats are supported
%   ANT - Advanced Neuro Technology, EEProbe (*.avr, *.eeg, *.cnt)
%   BCI2000 (*.dat)
%   Biosemi (*.bdf)
%   BrainVision (*.eeg, *.seg, *.dat, *.vhdr, *.vmrk)
%   CED - Cambridge Electronic Design (*.smr)
%   EGI - Electrical Geodesics, Inc. (*.egis, *.ave, *.gave, *.ses, *.raw, *.sbin, *.mff)
%   GTec (*.mat, *.hdf5)
%   Generic data formats (*.edf, *.gdf)
%   Megis/BESA (*.avr, *.swf, *.besa)
%   NeuroScan (*.eeg, *.cnt, *.avg)
%   Nexstim (*.nxe)
%   TMSi (*.Poly5)
%   Mega Neurone (directory)
%   Natus/Nicolet/Nervus (.e files)
%   Nihon Kohden (*.m00, *.EEG)
%   Bitalino OpenSignals (*.txt)
%   OpenBCI (*.txt)
%
% The following spike and LFP dataformats are supported
%   Neuralynx (*.ncs, *.nse, *.nts, *.nev, *.nrd, *.dma, *.log)
%   Plextor (*.nex, *.plx, *.ddt)
%   CED - Cambridge Electronic Design (*.smr)
%   MPI - Max Planck Institute (*.dap)
%   Neurosim  (neurosim_spikes, neurosim_signals, neurosim_ds)
%   Windaq (*.wdq)
%   NeuroOmega (*.mat transformed from *.mpx)
%   Neurodata Without Borders: Neurophysiology (*.nwb)
%
% The following NIRS dataformats are supported
%   Artinis - Artinis Medical Systems B.V. (*.oxy3, *.oxy4, *.oxyproj)
%   BUCN - Birkbeck college, London (*.txt)
%   SNIRF - Society for functional near-infrared spectroscopy (*.snirf)
%
% The following Eyetracker dataformats are supported
%   EyeLink - SR Research (*.asc)
%   SensoMotoric Instruments - (*.txt)
%   Tobii - (*.tsv)
%
% See also FT_READ_DATA, FT_READ_EVENT, FT_WRITE_DATA, FT_WRITE_EVENT,
% FT_CHANTYPE, FT_CHANUNIT

% Copyright (C) 2003-2021 Robert Oostenveld
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

% TODO channel renaming should be made a general option (see bham_bdf)

persistent cacheheader        % for caching the full header
persistent cachechunk         % for caching the res4 chunk when doing realtime analysis on the CTF scanner
persistent db_blob            % for fcdc_mysql

if isempty(db_blob)
  db_blob = false;
end

if iscell(filename)
  % use recursion to read the header from multiple files
  ft_warning('concatenating header from %d files', numel(filename));
  
  hdr = cell(size(filename));
  for i=1:numel(filename)
    hdr{i} = ft_read_header(filename{i}, varargin{:});
  end
  
  allhdr = cat(1, hdr{:});
  if numel(unique([allhdr.label]))==sum([allhdr.nChans])
    % each file has different channels, concatenate along the channel dimension
    for i=1:numel(filename)
      assert(isequal(hdr{i}.Fs, hdr{1}.Fs), 'sampling rates are not consistent over files');
      assert(isequal(hdr{i}.nSamples, hdr{1}.nSamples), 'number of samples is not consistent over files');
      assert(isequal(hdr{i}.nTrials, hdr{1}.nTrials), 'number of trials is not consistent over files');
    end
    combined          = hdr{1}; % copy the first header as the general one
    combined.label    = [allhdr.label];
    combined.chanunit = [allhdr.chanunit];
    combined.chantype = [allhdr.chantype];
    combined.nChans   = sum([allhdr.nChans]);
    combined.orig     = hdr;    % store the original header details of all files
  else
    % each file has the same channels, concatenate along the time dimension
    ntrl = nan(size(filename));
    nsmp = nan(size(filename));
    for i=1:numel(filename)
      assert(isequal(hdr{i}.Fs, hdr{1}.Fs), 'sampling rates are not consistent over files');
      assert(isequal(hdr{i}.label, hdr{1}.label), 'channels are not consistent over files');
      ntrl(i) = hdr{i}.nTrials;
      nsmp(i) = hdr{i}.nSamples;
    end
    % the subsequent code concatenates the files over time, i.e. each file has the same channels
    combined      = hdr{1}; % copy the first header as the general one
    combined.orig = hdr;    % store the original header details of all files
    if all(ntrl==1)
      % each file is a continuous recording
      combined.nTrials  = ntrl(1);
      combined.nSamples = sum(nsmp);
    elseif all(nsmp==nsmp(1))
      % each file holds segments of the same length
      combined.nTrials  = sum(ntrl);
      combined.nSamples = nsmp(1);
    else
      ft_error('cannot concatenate files');
    end
  end % concatenate over channels or over time
  % return the header of the concatenated datafiles
  hdr = combined;
  return
end

% get the options
headerformat   = ft_getopt(varargin, 'headerformat');
retry          = ft_getopt(varargin, 'retry', false);     % the default is not to retry reading the header
chanindx       = ft_getopt(varargin, 'chanindx');         % this is used for EDF with different sampling rates
coordsys       = ft_getopt(varargin, 'coordsys', 'head'); % this is used for ctf and neuromag_mne, it can be head or dewar
coilaccuracy   = ft_getopt(varargin, 'coilaccuracy');     % empty, or a number between 0-2
chantype       = ft_getopt(varargin, 'chantype', {});
password       = ft_getopt(varargin, 'password', struct([]));
readbids       = ft_getopt(varargin, 'readbids', 'ifmakessense');

% this should be a cell array
if ~iscell(chantype); chantype = {chantype}; end

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

if isempty(headerformat)
  % only do the autodetection if the format was not specified
  headerformat = ft_filetype(filename);
end

if iscell(headerformat)
  % this happens for datasets specified as cell-array for concatenation
  headerformat = headerformat{1};
end

if strcmp(headerformat, 'compressed')
  % the file is compressed, unzip on the fly
  inflated     = true;
  filename     = inflate_file(filename);
  headerformat = ft_filetype(filename);
else
  inflated     = false;
end

% for backward compatibility with https://github.com/fieldtrip/fieldtrip/issues/1585
if islogical(readbids)
  % it should be either yes/no/ifmakessense
  if readbids
    readbids = 'yes';
  else
    readbids = 'no';
  end
end

realtime = any(strcmp(headerformat, {'fcdc_buffer', 'ctf_shm', 'fcdc_mysql'}));

% The checkUniqueLabels flag is used for the realtime buffer in case
% it contains fMRI data. It prevents 1000000 voxel names to be checked
% for uniqueness. fMRI users will probably never use channel names
% for anything.

if realtime
  % skip the rest of the initial checks to increase the speed for realtime operation
  
  checkUniqueLabels = false;
  % the cache and fallback option should always be false for realtime processing
  cache    = false;
  fallback = false;
  
else
  % check whether the file or directory exists
  if  ~exist(filename, 'file')
    ft_error('file or directory ''%s'' does not exist', filename);
  end
  
  checkUniqueLabels = true;
  % get the rest of the options, this is skipped for realtime operation
  cache          = ft_getopt(varargin, 'cache');
  fallback       = ft_getopt(varargin, 'fallback');
  checkmaxfilter = ft_getopt(varargin, 'checkmaxfilter', true);
  
  if isempty(cache)
    if any(strcmp(headerformat, {'bci2000_dat', 'eyelink_asc', 'gtec_mat', 'gtec_hdf5', 'mega_neurone', 'nihonkohden_m00', 'smi_txt', 'biosig'}))
      cache = true;
    else
      cache = false;
    end
  end
  
  % ensure that the headerfile and datafile are defined, which are sometimes different than the name of the dataset
  [filename, headerfile, datafile] = dataset2files(filename, headerformat);
  if ~strcmp(filename, headerfile) && ~ft_filetype(filename, 'ctf_ds') && ~ft_filetype(filename, 'fcdc_buffer_offline') && ~ft_filetype(filename, 'fcdc_matbin')
    filename     = headerfile;                % this function should read the headerfile, not the dataset
    headerformat = ft_filetype(filename);     % update the filetype
  end
end % if skip initial check

% implement the caching in a data-format independent way
if cache && exist(headerfile, 'file') && ~isempty(cacheheader)
  % try to get the header from cache
  details = dir(headerfile);
  if isequal(details, cacheheader.details)
    % the header file has not been updated, fetch it from the cache
    % fprintf('got header from cache\n');
    hdr = rmfield(cacheheader, 'details');
    
    switch ft_filetype(datafile)
      case {'ctf_ds' 'ctf_meg4' 'ctf_old' 'read_ctf_res4'}
        % for realtime analysis end-of-file-chasing the res4 does not correctly
        % estimate the number of samples, so we compute it on the fly
        sz = 0;
        files = dir([filename '/*.*meg4']);
        for j=1:numel(files)
          sz = sz + files(j).bytes;
        end
        hdr.nTrials = floor((sz - 8) / (hdr.nChans*4) / hdr.nSamples);
    end
    
    return
  end % if the details correspond
end % if cache

% the support for head/dewar coordinates is still limited
if strcmp(coordsys, 'dewar') && ~any(strcmp(headerformat, {'fcdc_buffer', 'ctf_ds', 'ctf_meg4', 'ctf_res4', 'neuromag_fif', 'neuromag_mne'}))
  ft_error('dewar coordinates are not supported for %s', headerformat);
end

% deal with data that is organized according to BIDS
if strcmp(readbids, 'yes') || strcmp(readbids, 'ifmakessense')
  [p, f, x] = fileparts(filename);
  % check whether it is a BIDS dataset with json and tsv sidecar files
  % data in a BIDS tsv file (like physio and stim) will be explicitly dealt with in BIDS_TSV
  isbids = startsWith(f, 'sub-') && ~strcmp(x, '.tsv');
  if isbids
    % try to read the metadata from the BIDS sidecar files
    sidecar = bids_sidecar(filename);
    if ~isempty(sidecar)
      data_json = read_json(sidecar);
    end
    sidecar = bids_sidecar(filename, 'channels');
    if ~isempty(sidecar)
      channels_tsv = read_tsv(sidecar);
    end
    sidecar = bids_sidecar(filename, 'electrodes');
    if ~isempty(sidecar)
      electrodes_tsv = read_tsv(sidecar);
    end
    sidecar = bids_sidecar(filename, 'optodes');
    if ~isempty(sidecar)
      optodes_tsv = read_tsv(sidecar);
    end
  end
end

% start with an empty header
hdr = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the data with the low-level reading function
% please maintain this list in alphabetical order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch headerformat
  
  case '4d'
    orig            = read_4d_hdr(datafile);
    hdr.Fs          = orig.header_data.SampleFrequency;
    hdr.nChans      = orig.header_data.TotalChannels;
    hdr.nSamples    = orig.header_data.SlicesPerEpoch;
    hdr.nSamplesPre = round(orig.header_data.FirstLatency*orig.header_data.SampleFrequency);
    hdr.nTrials     = orig.header_data.TotalEpochs;
    %hdr.label       = {orig.channel_data(:).chan_label}';
    hdr.label       = orig.Channel;
    [hdr.grad, elec] = bti2grad(orig);
    if ~isempty(elec)
      hdr.elec = elec;
    end
    
    % remember original header details
    hdr.orig        = orig;
    
  case {'4d_pdf', '4d_m4d', '4d_xyz'}
    orig            = read_bti_m4d(filename);
    hdr.Fs          = orig.SampleFrequency;
    hdr.nChans      = orig.TotalChannels;
    hdr.nSamples    = orig.SlicesPerEpoch;
    hdr.nSamplesPre = round(orig.FirstLatency*orig.SampleFrequency);
    hdr.nTrials     = orig.TotalEpochs;
    hdr.label       = orig.ChannelOrder(:);
    [hdr.grad, elec] = bti2grad(orig);
    if ~isempty(elec)
      hdr.elec = elec;
    end
    
    % remember original header details
    hdr.orig        = orig;
    
  case 'AnyWave'
    orig = read_ahdf5_hdr(datafile);
    hdr.orig = orig;
    hdr.Fs = orig.channels(1).samplingRate;
    hdr.nChans = numel(orig.channels);
    hdr.nSamples = orig.numberOfSamples;
    hdr.nTrials = orig.numberOfBlocks;
    hdr.nSamplesPre = 0;
    hdr.label = orig.label;
    hdr.reference = orig.reference(:);
    hdr.chanunit = orig.unit(:);
    hdr.chantype = orig.type(:);
    
  case 'bci2000_dat'
    % this requires the load_bcidat mex file to be present on the path
    ft_hastoolbox('BCI2000', 1);
    % this is inefficient, since it reads the complete data
    [signal, states, parameters, total_samples] = load_bcidat(filename);
    % convert into a FieldTrip-like header
    hdr             = [];
    hdr.nChans      = size(signal,2);
    hdr.nSamples    = total_samples;
    hdr.nSamplesPre = 0;  % it is continuous
    hdr.nTrials     = 1;  % it is continuous
    hdr.Fs          = parameters.SamplingRate.NumericValue;
    % there are some differences in the fields that are present in the
    % *.dat files, probably due to different BCI2000 versions
    if isfield(parameters, 'ChannelNames') && isfield(parameters.ChannelNames, 'Value') && ~isempty(parameters.ChannelNames.Value)
      hdr.label       = parameters.ChannelNames.Value;
    elseif isfield(parameters, 'ChannelNames') && isfield(parameters.ChannelNames, 'Values') && ~isempty(parameters.ChannelNames.Values)
      hdr.label       = parameters.ChannelNames.Values;
    else
      % give this warning only once
      ft_warning('creating fake channel names');
      for i=1:hdr.nChans
        hdr.label{i} = sprintf('%d', i);
      end
    end
    
    % remember the original header details
    hdr.orig.parameters       = parameters;
    % also remember the complete data upon request
    if cache
      hdr.orig.signal         = signal;
      hdr.orig.states         = states;
      hdr.orig.total_samples  = total_samples;
    end
    
  case 'besa_besa'
    if isempty(chanindx)
      hdr = read_besa_besa(filename);
    else
      hdr = read_besa_besa(filename,[],1);
      if chanindx > hdr.orig.channel_info.orig_n_channels
        ft_error('FILEIO:InvalidChanIndx', 'selected channels are not present in the data');
      else
        hdr = read_besa_besa(filename,[],chanindx);
      end
    end
    
  case 'besa_avr'
    orig = read_besa_avr(filename);
    hdr.Fs          = 1000/orig.di;
    hdr.nChans      = size(orig.data,1);
    hdr.nSamples    = size(orig.data,2);
    hdr.nSamplesPre = -(hdr.Fs * orig.tsb/1000);   % convert from ms to samples
    hdr.nTrials     = 1;
    if isfield(orig, 'label') && iscell(orig.label)
      hdr.label = orig.label;
    elseif isfield(orig, 'label') && ischar(orig.label)
      hdr.label = tokenize(orig.label, ' ');
    else
      % give this warning only once
      ft_warning('creating fake channel names');
      for i=1:hdr.nChans
        hdr.label{i} = sprintf('%d', i);
      end
    end
    
  case 'besa_swf'
    orig = read_besa_swf(filename);
    hdr.Fs          = 1000/orig.di;
    hdr.nChans      = size(orig.data,1);
    hdr.nSamples    = size(orig.data,2);
    hdr.nSamplesPre = -(hdr.Fs * orig.tsb/1000);   % convert from ms to samples
    hdr.nTrials     = 1;
    hdr.label       = orig.label;
    
  case 'biosig'
    % this requires the biosig toolbox
    ft_hastoolbox('BIOSIG', 1);
    hdr = read_biosig_header(filename);
    
  case {'biosemi_bdf', 'bham_bdf'}
    hdr = read_biosemi_bdf(filename);
    if any(diff(hdr.orig.SampleRate))
      ft_error('channels with different sampling rate not supported');
    end
    
    if ~ft_senstype(hdr, 'ext1020')
      % assign the channel type and units for the known channels
      hdr.chantype = repmat({'unknown'}, size(hdr.label));
      hdr.chanunit = repmat({'unknown'}, size(hdr.label));
      chan = ~cellfun(@isempty, regexp(hdr.label, '^[A-D]\d*$'));
      hdr.chantype(chan) = {'eeg'};
      hdr.chanunit(chan) = {'uV'};
    end
    
    if ft_filetype(filename, 'bham_bdf')
      % TODO channel renaming should be made a general option
      % this is for the Biosemi system used at the University of Birmingham
      labelold = { 'A1' 'A2' 'A3' 'A4' 'A5' 'A6' 'A7' 'A8' 'A9' 'A10' 'A11' 'A12' 'A13' 'A14' 'A15' 'A16' 'A17' 'A18' 'A19' 'A20' 'A21' 'A22' 'A23' 'A24' 'A25' 'A26' 'A27' 'A28' 'A29' 'A30' 'A31' 'A32' 'B1' 'B2' 'B3' 'B4' 'B5' 'B6' 'B7' 'B8' 'B9' 'B10' 'B11' 'B12' 'B13' 'B14' 'B15' 'B16' 'B17' 'B18' 'B19' 'B20' 'B21' 'B22' 'B23' 'B24' 'B25' 'B26' 'B27' 'B28' 'B29' 'B30' 'B31' 'B32' 'C1' 'C2' 'C3' 'C4' 'C5' 'C6' 'C7' 'C8' 'C9' 'C10' 'C11' 'C12' 'C13' 'C14' 'C15' 'C16' 'C17' 'C18' 'C19' 'C20' 'C21' 'C22' 'C23' 'C24' 'C25' 'C26' 'C27' 'C28' 'C29' 'C30' 'C31' 'C32' 'D1' 'D2' 'D3' 'D4' 'D5' 'D6' 'D7' 'D8' 'D9' 'D10' 'D11' 'D12' 'D13' 'D14' 'D15' 'D16' 'D17' 'D18' 'D19' 'D20' 'D21' 'D22' 'D23' 'D24' 'D25' 'D26' 'D27' 'D28' 'D29' 'D30' 'D31' 'D32' 'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8' 'Status'};
      labelnew = { 'P9' 'PPO9h' 'PO7' 'PPO5h' 'PPO3h' 'PO5h' 'POO9h' 'PO9' 'I1' 'OI1h' 'O1' 'POO1' 'PO3h' 'PPO1h' 'PPO2h' 'POz' 'Oz' 'Iz' 'I2' 'OI2h' 'O2' 'POO2' 'PO4h' 'PPO4h' 'PO6h' 'POO10h' 'PO10' 'PO8' 'PPO6h' 'PPO10h' 'P10' 'P8' 'TPP9h' 'TP7' 'TTP7h' 'CP5' 'TPP7h' 'P7' 'P5' 'CPP5h' 'CCP5h' 'CP3' 'P3' 'CPP3h' 'CCP3h' 'CP1' 'P1' 'Pz' 'CPP1h' 'CPz' 'CPP2h' 'P2' 'CPP4h' 'CP2' 'CCP4h' 'CP4' 'P4' 'P6' 'CPP6h' 'CCP6h' 'CP6' 'TPP8h' 'TP8' 'TPP10h' 'T7' 'FTT7h' 'FT7' 'FC5' 'FCC5h' 'C5' 'C3' 'FCC3h' 'FC3' 'FC1' 'C1' 'CCP1h' 'Cz' 'FCC1h' 'FCz' 'FFC1h' 'Fz' 'FFC2h' 'FC2' 'FCC2h' 'CCP2h' 'C2' 'C4' 'FCC4h' 'FC4' 'FC6' 'FCC6h' 'C6' 'TTP8h' 'T8' 'FTT8h' 'FT8' 'FT9' 'FFT9h' 'F7' 'FFT7h' 'FFC5h' 'F5' 'AFF7h' 'AF7' 'AF5h' 'AFF5h' 'F3' 'FFC3h' 'F1' 'AF3h' 'Fp1' 'Fpz' 'Fp2' 'AFz' 'AF4h' 'F2' 'FFC4h' 'F4' 'AFF6h' 'AF6h' 'AF8' 'AFF8h' 'F6' 'FFC6h' 'FFT8h' 'F8' 'FFT10h' 'FT10'};
      % rename the channel labels
      for i=1:length(labelnew)
        chan = strcmp(labelold(i), hdr.label);
        hdr.label(chan) = labelnew(chan);
      end
    end
    
  case {'biosemi_old'}
    % this uses the openbdf and readbdf functions that were copied from EEGLAB
    orig = openbdf(filename);
    if any(orig.Head.SampleRate~=orig.Head.SampleRate(1))
      ft_error('channels with different sampling rate not supported');
    end
    hdr.Fs          = orig.Head.SampleRate(1);
    hdr.nChans      = orig.Head.NS;
    hdr.label       = cellstr(orig.Head.Label);
    % it is continuous data, therefore append all records in one trial
    hdr.nSamples    = orig.Head.NRec * orig.Head.Dur * orig.Head.SampleRate(1);
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1;
    hdr.orig        = orig;
    % close the file between separate read operations
    fclose(orig.Head.FILE.FID);
    
  case 'blackrock_nev'
    % read header for nsX file associated with NEV file
    % use ft_read_event to read event information in .nev file
    
    ft_hastoolbox('NPMK', 1);
    % ensure that the filename contains a full path specification,
    % otherwise the low-level function fails
    [p,n] = fileparts(filename);
    if isempty(p)
      filename = which(filename);
      [p,n] = fileparts(filename);
    end
    
    NEV = openNEV(filename,'noread','nosave');
    
    %searching for associated nsX file in same folder
    files=dir(strcat(fullfile(p,n),'.ns*'));
    if isempty(files)
      ft_error('no .ns* file associated to %s in %s',n,p);
    end
    
    %searching for nsX file with same sampling freq that NEV
    for i=1:numel(files)
      nsX_hdr = ft_read_header(fullfile(p,files(i).name),'chantype',chantype);
      if nsX_hdr.Fs == NEV.MetaTags.SampleRes
        hdr = nsX_hdr;
        break
      end
    end
    
    if isempty(hdr)
      ft_error('no .ns* file with same sampling frequency as %s (%i)',n,NEV.MetaTags.SampleRes);
    end
    
  case 'blackrock_nsx'
    ft_hastoolbox('NPMK', 1);
    % ensure that the filename contains a full path specification,
    % otherwise the low-level function fails
    p = fileparts(filename);
    if isempty(p)
      filename = which(filename);
    end
    
    orig = openNSx(filename, 'noread');
    channelstype=regexp({orig.ElectrodesInfo.Label},'[A-Za-z]+','match','once');
    chaninfo=table({orig.ElectrodesInfo.ElectrodeID}',...
      transpose(deblank({orig.ElectrodesInfo.Label})),[channelstype]',...
      {orig.ElectrodesInfo.ConnectorBank}',{orig.ElectrodesInfo.ConnectorPin}',...
      transpose(deblank({orig.ElectrodesInfo.AnalogUnits})),...
      'VariableNames',{'id' 'label' 'chantype' 'bank' 'pin' 'unit'});
    
    if isempty(chantype)
      chantype = unique(channelstype,'stable');
    end
    
    %selecting channel according to chantype
    orig_label=deblank({orig.ElectrodesInfo.Label});
    orig_unit=deblank({orig.ElectrodesInfo.AnalogUnits});
    channels={}; channelstype={}; channelsunit={}; skipfactor=[];
    for c=1:length(chantype)
      chantype_split=strsplit(chantype{c},':');
      if numel(chantype_split) == 2
        chantype{c}=chantype_split{1};
        skipfactor=[skipfactor,str2double(chantype_split{2})];
      elseif numel(chantype_split) > 2
        ft_error('Use : to specify skipfactor, e.g. analog:10')
      end
      chan_sel=~cellfun(@isempty,regexp(orig_label,chantype{c}));
      if sum(chan_sel)==0
        if ~strcmp(chantype{c},'chaninfo')
          ft_error('unknown chantype %s, available channels are %s',chantype{c},strjoin(orig_label));
        end
      else
        channels=[channels, orig_label(chan_sel)];
        channelsunit=[channelsunit, orig_unit(chan_sel)];
        channelstype=[channelstype, repmat(chantype(c), 1, sum(chan_sel))];
      end
    end
    
    skipfactor=unique(skipfactor);
    if isempty(skipfactor)
      skipfactor=1;
    elseif length(skipfactor)>1
      ft_error('inconsistent skip factors across channels');
    end
    
    %If no channel selected issue error specifying available chantypes
    if isempty(channels)
      ft_error('No channel selected. Availabe chantypes are: %s',strjoin(unique(chaninfo.chantype)));
    end
    
    hdr.Fs          = orig.MetaTags.SamplingFreq/skipfactor;
    hdr.nChans      = length(channels);
    hdr.nSamples    = floor(orig.MetaTags.DataPoints/skipfactor);
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1; %?
    hdr.label       = deblank(channels)';
    hdr.chantype    = channelstype;
    hdr.chanunit    = channelsunit;
    hdr.orig        = orig;
    hdr.orig.chaninfo = chaninfo;
    hdr.orig.skipfactor = skipfactor;
    
  case {'brainvision_vhdr', 'brainvision_seg', 'brainvision_eeg', 'brainvision_dat'}
    orig = read_brainvision_vhdr(filename);
    hdr.Fs          = orig.Fs;
    hdr.nChans      = orig.NumberOfChannels;
    hdr.label       = orig.label;
    hdr.nSamples    = orig.nSamples;
    hdr.nSamplesPre = orig.nSamplesPre;
    hdr.nTrials     = orig.nTrials;
    hdr.orig        = orig;
    % assign the channel type and units for the known channels
    hdr.chantype = repmat({'eeg'}, size(hdr.label));
    hdr.chanunit = repmat({'uV'},  size(hdr.label));
    
  case 'bucn_nirs'
    orig = read_bucn_nirshdr(filename);
    hdr  = rmfield(orig, 'time');
    hdr.orig = orig;
    
  case 'ced_son'
    % check that the required low-level toolbox is available
    ft_hastoolbox('neuroshare', 1);
    % use the reading function supplied by Gijs van Elswijk
    orig = read_ced_son(filename,'readevents','no','readdata','no');
    orig = orig.header;
    % In Spike2, channels can have different sampling rates, units, length
    % etc. etc. Here, channels need to have to same properties.
    if length(unique([orig.samplerate]))>1,
      ft_error('channels with different sampling rates are not supported');
    else
      hdr.Fs   = orig(1).samplerate;
    end
    hdr.nChans = length(orig);
    % nsamples of the channel with least samples
    hdr.nSamples    = min([orig.nsamples]);
    hdr.nSamplesPre = 0;
    % only continuous data supported
    if sum(strcmpi({orig.mode},'continuous')) < hdr.nChans,
      ft_error('not all channels contain continuous data');
    else
      hdr.nTrials = 1;
    end
    hdr.label = {orig.label};
    
  case  'combined_ds'
    hdr = read_combined_ds(filename);
    
  case {'ctf_ds', 'ctf_meg4', 'ctf_res4'}
    % check the presence of the required low-level toolbox
    ft_hastoolbox('ctf', 1);
    orig             = readCTFds(filename);
    if isempty(orig)
      % this is to deal with data from the 64 channel system and the error
      % readCTFds: .meg4 file header=MEG4CPT   Valid header options:  MEG41CP  MEG42CP
      ft_error('could not read CTF with this implementation, please try again with the ''ctf_old'' file format');
    end
    hdr.Fs           = orig.res4.sample_rate;
    hdr.nChans       = orig.res4.no_channels;
    hdr.nSamples     = orig.res4.no_samples;
    hdr.nSamplesPre  = orig.res4.preTrigPts;
    hdr.nTrials      = orig.res4.no_trials;
    hdr.label        = cellstr(orig.res4.chanNames);
    for i=1:numel(hdr.label)
      % remove the site-specific numbers from each channel name, e.g. 'MZC01-1706' becomes 'MZC01'
      hdr.label{i} = strtok(hdr.label{i}, '-');
    end
    % read the balance coefficients, these are used to compute the synthetic gradients
    coeftype = cellstr(char(orig.res4.scrr(:).coefType));
    try
      [alphaMEG,MEGlist,Refindex] = getCTFBalanceCoefs(orig,'NONE', 'T');
      orig.BalanceCoefs.none.alphaMEG  = alphaMEG;
      orig.BalanceCoefs.none.MEGlist   = MEGlist;
      orig.BalanceCoefs.none.Refindex  = Refindex;
    catch
      ft_warning('cannot read balancing coefficients for NONE');
    end
    if any(~cellfun(@isempty,strfind(coeftype, 'G1BR')))
      try
        [alphaMEG,MEGlist,Refindex] = getCTFBalanceCoefs(orig,'G1BR', 'T');
        orig.BalanceCoefs.G1BR.alphaMEG  = alphaMEG;
        orig.BalanceCoefs.G1BR.MEGlist   = MEGlist;
        orig.BalanceCoefs.G1BR.Refindex  = Refindex;
      catch
        ft_warning('cannot read balancing coefficients for G1BR');
      end
    end
    if any(~cellfun(@isempty,strfind(coeftype, 'G2BR')))
      try
        [alphaMEG,MEGlist,Refindex] = getCTFBalanceCoefs(orig, 'G2BR', 'T');
        orig.BalanceCoefs.G2BR.alphaMEG  = alphaMEG;
        orig.BalanceCoefs.G2BR.MEGlist   = MEGlist;
        orig.BalanceCoefs.G2BR.Refindex  = Refindex;
      catch
        ft_warning('cannot read balancing coefficients for G2BR');
      end
    end
    if any(~cellfun(@isempty,strfind(coeftype, 'G3BR')))
      try
        [alphaMEG,MEGlist,Refindex] = getCTFBalanceCoefs(orig, 'G3BR', 'T');
        orig.BalanceCoefs.G3BR.alphaMEG  = alphaMEG;
        orig.BalanceCoefs.G3BR.MEGlist   = MEGlist;
        orig.BalanceCoefs.G3BR.Refindex  = Refindex;
      catch
        ft_warning('cannot read balancing coefficients for G3BR');
      end
    end
    if any(~cellfun(@isempty,strfind(coeftype, 'G3AR')))
      try
        [alphaMEG,MEGlist,Refindex] = getCTFBalanceCoefs(orig, 'G3AR', 'T');
        orig.BalanceCoefs.G3AR.alphaMEG  = alphaMEG;
        orig.BalanceCoefs.G3AR.MEGlist   = MEGlist;
        orig.BalanceCoefs.G3AR.Refindex  = Refindex;
      catch
        % May not want a warning here if these are not commonly used.
        % Already get a (fprintf) warning from getCTFBalanceCoefs.m
        % ft_warning('cannot read balancing coefficients for G3AR');
      end
    end
    % add a gradiometer structure for forward and inverse modelling
    try
      [grad, elec] = ctf2grad(orig, strcmp(coordsys, 'dewar'), coilaccuracy);
      if ~isempty(grad)
        hdr.grad = grad;
      end
      if ~isempty(elec)
        hdr.elec = elec;
      end
    catch
      % this fails if the res4 file is not correctly closed, e.g. during realtime processing
      tmp = lasterror;
      disp(tmp.message);
      ft_warning('could not construct gradiometer definition from the header');
    end
    
    % for realtime analysis EOF chasing the res4 does not correctly
    % estimate the number of samples, so we compute it on the fly from the
    % meg4 file sizes.
    sz = 0;
    files = dir([filename '/*.*meg4']);
    for j=1:numel(files)
      sz = sz + files(j).bytes;
    end
    hdr.nTrials = floor((sz - 8) / (hdr.nChans*4) / hdr.nSamples);
    
    % add the original header details
    hdr.orig = orig;
    
  case {'ctf_old', 'read_ctf_res4'}
    % read it using the open-source MATLAB code that originates from CTF and that was modified by the FCDC
    orig             = read_ctf_res4(headerfile);
    hdr.Fs           = orig.Fs;
    hdr.nChans       = orig.nChans;
    hdr.nSamples     = orig.nSamples;
    hdr.nSamplesPre  = orig.nSamplesPre;
    hdr.nTrials      = orig.nTrials;
    hdr.label        = orig.label;
    % add a gradiometer structure for forward and inverse modelling
    try
      hdr.grad = ctf2grad(orig);
    catch
      % this fails if the res4 file is not correctly closed, e.g. during realtime processing
      tmp = lasterror;
      disp(tmp.message);
      ft_warning('could not construct gradiometer definition from the header');
    end
    % add the original header details
    hdr.orig = orig;
    
  case 'ctf_read_res4'
    % check that the required low-level toolbos ix available
    ft_hastoolbox('eegsf', 1);
    % read it using the CTF importer from the NIH and Daren Weber
    orig = ctf_read_res4(fileparts(headerfile), 0);
    % convert the header into a structure that FieldTrip understands
    hdr              = [];
    hdr.Fs           = orig.setup.sample_rate;
    hdr.nChans       = length(orig.sensor.info);
    hdr.nSamples     = orig.setup.number_samples;
    hdr.nSamplesPre  = orig.setup.pretrigger_samples;
    hdr.nTrials      = orig.setup.number_trials;
    for i=1:length(orig.sensor.info)
      hdr.label{i}   = orig.sensor.info(i).label;
    end
    hdr.label        = hdr.label(:);
    % add a gradiometer structure for forward and inverse modelling
    try
      hdr.grad = ctf2grad(orig);
    catch
      % this fails if the res4 file is not correctly closed, e.g. during realtime processing
      tmp = lasterror;
      disp(tmp.message);
      ft_warning('could not construct gradiometer definition from the header');
    end
    % add the original header details
    hdr.orig = orig;
    
  case 'ctf_shm'
    % read the header information from shared memory
    hdr = read_shm_header(filename);
    
    
  case {'curry_dat', 'curry_cdt'}
    orig            = load_curry_data_file(filename);
    hdr             = [];
    hdr.Fs          = orig.fFrequency;
    hdr.nChans      = orig.nChannels;
    hdr.nSamples    = orig.nSamples;
    hdr.nSamplesPre = sum(orig.time<0);
    hdr.nTrials     = orig.nTrials;
    hdr.label       = orig.labels(:);
    hdr.orig        = orig;
    
  case 'dataq_wdq'
    orig            = read_wdq_header(filename);
    hdr             = [];
    hdr.Fs          = orig.fsample;
    hdr.nChans      = orig.nchan;
    hdr.nSamples    = orig.nbytesdat/(2*hdr.nChans);
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1;
    for k = 1:hdr.nChans
      if isfield(orig.chanhdr(k), 'annot') && ~isempty(orig.chanhdr(k).annot)
        hdr.label{k,1} = orig.chanhdr(k).annot;
      else
        hdr.label{k,1} = orig.chanhdr(k).label;
      end
    end
    
    % add the original header details
    hdr.orig  = orig;
    
  case {'deymed_ini' 'deymed_dat'}
    % the header is stored in a *.ini file
    orig            = read_deymed_ini(headerfile);
    hdr             = [];
    hdr.Fs          = orig.Fs;
    hdr.nChans      = orig.nChans;
    hdr.nSamples    = orig.nSamples;
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1;
    hdr.label       = orig.label(:);
    hdr.orig        = orig; % remember the original details
    
  case 'edf'
    % this reader is largely similar to the bdf reader
    if isempty(chanindx)
      hdr = read_edf(filename);
    else
      hdr = read_edf(filename,[],1);
      if chanindx > hdr.orig.NS
        ft_error('FILEIO:InvalidChanIndx', 'selected channels are not present in the data');
      else
        hdr = read_edf(filename,[],chanindx);
      end
    end
    
  case 'eep_avr'
    % check that the required low-level toolbox is available
    ft_hastoolbox('eeprobe', 1);
    % read the whole average and keep only header info (it is a bit silly, but the easiest to do here)
    hdr = read_eep_avr(filename);
    hdr.Fs          = hdr.rate;
    hdr.nChans      = size(hdr.data,1);
    hdr.nSamples    = size(hdr.data,2);
    hdr.nSamplesPre = hdr.xmin*hdr.rate/1000;
    hdr.nTrials     = 1;        % it can always be interpreted as continuous data
    % remove the data and variance if present
    hdr = removefields(hdr, {'data', 'variance'});
    
  case 'eep_cnt'
    % check that the required low-level toolbox is available
    ft_hastoolbox('eeprobe', 1);
    % read the first sample from the continuous data, this will also return the header
    orig = read_eep_cnt(filename, 1, 1);
    hdr.Fs          = orig.rate;
    hdr.nSamples    = orig.nsample;
    hdr.nSamplesPre = 0;
    hdr.label       = orig.label;
    hdr.nChans      = orig.nchan;
    hdr.nTrials     = 1;        % it can always be interpreted as continuous data
    hdr.orig        = orig;     % remember the original details
    
    
  case 'eeglab_set'
    hdr = read_eeglabheader(filename);
    
  case 'eeglab_erp'
    hdr = read_erplabheader(filename);
    
  case 'emotiv_mat'
    % This is a MATLAB *.mat file that is created using the Emotiv MATLAB
    % example code. It contains a 25xNsamples matrix and some other stuff.
    orig = load(filename);
    hdr.Fs          = 128;
    hdr.nChans      = 25;
    hdr.nSamples    = size(orig.data_eeg,1);
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1;
    hdr.label       = {'ED_COUNTER','ED_INTERPOLATED','ED_RAW_CQ','ED_AF3','ED_F7','ED_F3','ED_FC5','ED_T7','ED_P7','ED_O1','ED_O2','ED_P8','ED_T8','ED_FC6','ED_F4','ED_F8','ED_AF4','ED_GYROX','ED_GYROY','ED_TIMESTAMP','ED_ES_TIMESTAMP','ED_FUNC_ID','ED_FUNC_VALUE','ED_MARKER','ED_SYNC_SIGNAL'};
    % store the complete information in hdr.orig
    % ft_read_data and ft_read_event will get it from there
    hdr.orig        = orig;
    
  case 'eyelink_asc'
    asc = read_eyelink_asc(filename);
    hdr.nChans              = size(asc.dat,1);
    hdr.nSamples            = size(asc.dat,2);
    hdr.nSamplesPre         = 0;
    hdr.nTrials             = 1;
    hdr.FirstTimeStamp      = asc.dat(1,1);
    hdr.TimeStampPerSample  = mean(diff(asc.dat(1,:)));
    hdr.Fs                  = 1000/hdr.TimeStampPerSample;  % these timestamps are in miliseconds
    % give this warning only once
    ft_warning('creating fake channel names');
    for i=1:hdr.nChans
      hdr.label{i} = sprintf('%d', i);
    end
    
    % remember the original header details
    hdr.orig.header = asc.header;
    % remember all header and data details upon request
    if cache
      hdr.orig = asc;
    end
    
  case  'spmeeg_mat'
    hdr = read_spmeeg_header(filename);
    
  case  'ced_spike6mat'
    hdr = read_spike6mat_header(filename);
    
  case 'egi_egia'
    [fhdr,chdr,ename,cnames,fcom,ftext] = read_egis_header(filename);
    [p, f, x]       = fileparts(filename);
    
    if any(chdr(:,4)-chdr(1,4))
      ft_error('Sample rate not the same for all cells.');
    end
    
    hdr.Fs          = chdr(1,4); %making assumption that sample rate is same for all cells
    hdr.nChans      = fhdr(19);
    for i = 1:hdr.nChans
      % this should be consistent with ft_senslabel
      hdr.label{i,1}  = ['E' num2str(i)];
    end
    %since NetStation does not properly set the fhdr(11) field, use the number of subjects from the chdr instead
    hdr.nTrials     = chdr(1,2)*fhdr(18); %number of trials is numSubjects * numCells
    hdr.nSamplesPre = ceil(fhdr(14)/(1000/hdr.Fs));
    
    if any(chdr(:,3)-chdr(1,3))
      ft_error('Number of samples not the same for all cells.');
    end
    
    hdr.nSamples    = chdr(1,3); %making assumption that number of samples is same for all cells
    
    % remember the original header details
    hdr.orig.fhdr   = fhdr;
    hdr.orig.chdr   = chdr;
    hdr.orig.ename  = ename;
    hdr.orig.cnames = cnames;
    hdr.orig.fcom   = fcom;
    hdr.orig.ftext  = ftext;
    
  case 'egi_egis'
    [fhdr,chdr,ename,cnames,fcom,ftext] = read_egis_header(filename);
    [p, f, x]       = fileparts(filename);
    
    if any(chdr(:,4)-chdr(1,4))
      ft_error('Sample rate not the same for all cells.');
    end
    
    hdr.Fs          = chdr(1,4); %making assumption that sample rate is same for all cells
    hdr.nChans      = fhdr(19);
    for i = 1:hdr.nChans
      % this should be consistent with ft_senslabel
      hdr.label{i,1}  = ['E' num2str(i)];
    end
    hdr.nTrials     = sum(chdr(:,2));
    hdr.nSamplesPre = ceil(fhdr(14)/(1000/hdr.Fs));
    % assuming that a utility was used to insert the correct baseline
    % duration into the header since it is normally absent. This slot is
    % actually allocated to the age of the subject, although NetStation
    % does not use it when generating an EGIS session file.
    
    if any(chdr(:,3)-chdr(1,3))
      ft_error('Number of samples not the same for all cells.');
    end
    
    hdr.nSamples    = chdr(1,3); %making assumption that number of samples is same for all cells
    
    % remember the original header details
    hdr.orig.fhdr   = fhdr;
    hdr.orig.chdr   = chdr;
    hdr.orig.ename  = ename;
    hdr.orig.cnames = cnames;
    hdr.orig.fcom   = fcom;
    hdr.orig.ftext  = ftext;
    
  case 'egi_sbin'
    [header_array, CateNames, CatLengths, preBaseline] = read_sbin_header(filename);
    [p, f, x]       = fileparts(filename);
    
    hdr.Fs          = header_array(9);
    hdr.nChans      = header_array(10);
    for i = 1:hdr.nChans
      % this should be consistent with ft_senslabel
      hdr.label{i,1}  = ['E' num2str(i)];
    end
    hdr.nTrials     = header_array(15);
    hdr.nSamplesPre = preBaseline;
    
    hdr.nSamples    = header_array(16); % making assumption that number of samples is same for all cells
    
    % remember the original header details
    hdr.orig.header_array   = header_array;
    hdr.orig.CateNames   = CateNames;
    hdr.orig.CatLengths  = CatLengths;
    
  case 'egi_mff_v1'
    % The following represents the code that was written by Ingrid, Robert
    % and Giovanni to get started with the EGI mff dataset format. It might
    % not support all details of the file formats.
    %
    % An alternative implementation has been provided by EGI, this is
    % released as fieldtrip/external/egi_mff and referred further down in
    % this function as 'egi_mff_v2'.
    %
    % An more recent implementation has been provided by EGI and Arno Delorme, this
    % is released as https://github.com/arnodelorme/mffmatlabio and referred further
    % down in this function as 'egi_mff_v3'.
    
    if ~usejava('jvm')
      ft_error('the xml2struct requires MATLAB to be running with the Java virtual machine (JVM)');
      % an alternative implementation which does not require the JVM but runs much slower is
      % available from http://www.mathworks.com/matlabcentral/fileexchange/6268-xml4mat-v2-0
    end
    
    % get header info from .bin files
    binfiles = dir(fullfile(filename, 'signal*.bin'));
    if isempty(binfiles)
      ft_error('could not find any signal.bin in mff directory')
    end
    
    orig = [];
    for iSig = 1:length(binfiles)
      signalname = binfiles(iSig).name;
      fullsignalname = fullfile(filename, signalname);
      orig.signal(iSig).blockhdr = read_mff_bin(fullsignalname);
    end
    
    % get hdr info from xml files
    ws = ft_warning('off', 'MATLAB:REGEXP:deprecated'); % due to some small code xml2struct
    xmlfiles = dir( fullfile(filename, '*.xml'));
    disp('reading xml files to obtain header info...')
    for i = 1:numel(xmlfiles)
      if strcmpi(xmlfiles(i).name(1:2), '._') % Mac sometimes creates this useless files, don't use them
      elseif strcmpi(xmlfiles(i).name(1:6), 'Events') % don't read in events here, can take a lot of time, and we can do that in ft_read_event
      else
        fieldname     = xmlfiles(i).name(1:end-4);
        filename_xml  = fullfile(filename, xmlfiles(i).name);
        orig.xml.(fieldname) = xml2struct(filename_xml);
      end
    end
    ft_warning(ws); % revert the warning state
    
    % epochs.xml seems the most common version, but epoch.xml might also
    % occur, so use only one name
    if isfield(orig.xml, 'epoch')
      orig.xml.epochs = orig.xml.epoch;
      orig.xml = rmfield(orig.xml, 'epoch');
    end
    
    % make hdr according to FieldTrip rules
    hdr = [];
    Fs = zeros(length(orig.signal),1);
    nChans = zeros(length(orig.signal),1);
    nSamples = zeros(length(orig.signal),1);
    
    for iSig = 1:length(orig.signal)
      Fs(iSig)      = orig.signal(iSig).blockhdr(1).fsample(1);
      nChans(iSig)  = orig.signal(iSig).blockhdr(1).nsignals;
      % the number of samples per block can be different
      nSamples_Block = zeros(length(orig.signal(iSig).blockhdr),1);
      for iBlock  = 1:length(orig.signal(iSig).blockhdr)
        nSamples_Block(iBlock) = orig.signal(iSig).blockhdr(iBlock).nsamples(1);
      end
      nSamples(iSig) = sum(nSamples_Block);
    end
    
    if length(unique(Fs)) > 1 || length(unique(nSamples)) > 1
      ft_error('Fs and nSamples should be the same in all signals')
    end
    
    hdr.Fs          = Fs(1);
    hdr.nChans      = sum(nChans);
    hdr.nSamplesPre = 0;
    hdr.nSamples    = nSamples(1);
    hdr.nTrials     = 1;
    
    % get channel labels for signal 1 (main net), otherwise create them
    if isfield(orig.xml, 'sensorLayout') % assuming that signal1 is hdEEG sensornet, and channels are in xml file sensorLayout
      for iSens = 1:numel(orig.xml.sensorLayout.sensors)
        if ~isempty(orig.xml.sensorLayout.sensors(iSens).sensor.name) && ~(isstruct(orig.xml.sensorLayout.sensors(iSens).sensor.name) && numel(fieldnames(orig.xml.sensorLayout.sensors(iSens).sensor.name))==0)
          %only get name when channel is EEG (type 0), or REF (type 1),
          %rest are non interesting channels like place holders and COM and should not be added.
          if strcmp(orig.xml.sensorLayout.sensors(iSens).sensor.type, '0') || strcmp(orig.xml.sensorLayout.sensors(iSens).sensor.type, '1')
            % get the sensor name from the datafile
            hdr.label{iSens} = orig.xml.sensorLayout.sensors(iSens).sensor.name;
          end
        elseif strcmp(orig.xml.sensorLayout.sensors(iSens).sensor.type, '0') % EEG chan
          % this should be consistent with ft_senslabel
          hdr.label{iSens} = ['E' num2str(orig.xml.sensorLayout.sensors(iSens).sensor.number)];
        elseif strcmp(orig.xml.sensorLayout.sensors(iSens).sensor.type, '1') % REF chan
          % ingnie: I now choose REF as name for REF channel since our discussion see bug 1407. Arbitrary choice...
          hdr.label{iSens} = ['REF' num2str(iSens)];
        else
          % non interesting channels like place holders and COM
        end
      end
      % check if the amount of lables corresponds with nChannels in signal 1
      if length(hdr.label) == nChans(1)
        % good
      elseif length(hdr.label) > orig.signal(1).blockhdr(1).nsignals
        ft_warning('found more lables in xml.sensorLayout than channels in signal 1, thus can not use info in sensorLayout, creating labels on the fly')
        for iSens = 1:orig.signal(1).blockhdr(1).nsignals
          % this should be consistent with ft_senslabel
          hdr.label{iSens} = ['E' num2str(iSens)];
        end
      else
        ft_warning('found less lables in xml.sensorLayout than channels in signal 1, thus can not use info in sensorLayout, creating labels on the fly')
        for iSens = 1:orig.signal(1).blockhdr(1).nsignals
          % this should be consistent with ft_senslabel
          hdr.label{iSens} = ['E' num2str(iSens)];
        end
      end
      % get lables for other signals
      if length(orig.signal) == 2
        if isfield(orig.xml, 'pnsSet') % signal2 is PIB box, and lables are in xml file pnsSet
          nbEEGchan = length(hdr.label);
          for iSens = 1:numel(orig.xml.pnsSet.sensors)
            hdr.label{nbEEGchan+iSens} = num2str(orig.xml.pnsSet.sensors(iSens).sensor.name);
          end
          if length(hdr.label) == orig.signal(1).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
            % good
          elseif length(hdr.label) < orig.signal(1).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
            ft_warning('found less lables in xml.pnsSet than channels in signal 2, labeling with s2_unknownN instead')
            for iSens = length(hdr.label)+1 : orig.signal(1).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
              hdr.label{iSens} = ['s2_unknown', num2str(iSens)];
            end
          else
            ft_warning('found more lables in xml.pnsSet than channels in signal 2, thus can not use info in pnsSet, and labeling with s2_eN instead')
            for iSens = orig.signal(1).blockhdr(1).nsignals+1 : orig.signal(1).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
              hdr.label{iSens} = ['s2_E', num2str(iSens)];
            end
          end
        else % signal2 is not PIBbox
          ft_warning('creating channel labels for signal 2 on the fly')
          for iSens = 1:orig.signal(2).blockhdr(1).nsignals
            hdr.label{end+1} = ['s2_E', num2str(iSens)];
          end
        end
      elseif length(orig.signal) > 2
        % loop over signals and label channels accordingly
        ft_warning('creating channel labels for signal 2 to signal N on the fly')
        for iSig = 2:length(orig.signal)
          for iSens = 1:orig.signal(iSig).blockhdr(1).nsignals
            if iSig == 1 && iSens == 1
              hdr.label{1} = ['s',num2str(iSig),'_E', num2str(iSens)];
            else
              hdr.label{end+1} = ['s',num2str(iSig),'_E', num2str(iSens)];
            end
          end
        end
      end
    else % no xml.sensorLayout present
      ft_warning('no sensorLayout found in xml files, creating channel labels on the fly')
      for iSig = 1:length(orig.signal)
        for iSens = 1:orig.signal(iSig).blockhdr(1).nsignals
          if iSig == 1 && iSens == 1
            hdr.label{1} = ['s',num2str(iSig),'_E', num2str(iSens)];
          else
            hdr.label{end+1} = ['s',num2str(iSig),'_E', num2str(iSens)];
          end
        end
      end
    end
    
    % check if multiple epochs are present
    if isfield(orig.xml,'epochs')
      % add info to header about which sample correspond to which epochs, becasue this is quite hard for user to get...
      epochdef = zeros(length(orig.xml.epochs),3);
      for iEpoch = 1:length(orig.xml.epochs)
        if iEpoch == 1
          epochdef(iEpoch,1) = round(str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./(1000000./hdr.Fs))+1;
          epochdef(iEpoch,2) = round(str2double(orig.xml.epochs(iEpoch).epoch.endTime  )./(1000000./hdr.Fs));
          epochdef(iEpoch,3) = round(str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./(1000000./hdr.Fs)); % offset corresponds to timing
        else
          NbSampEpoch = round(str2double(orig.xml.epochs(iEpoch).epoch.endTime)./(1000000./hdr.Fs) - str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./(1000000./hdr.Fs));
          epochdef(iEpoch,1) = epochdef(iEpoch-1,2) + 1;
          epochdef(iEpoch,2) = epochdef(iEpoch-1,2) + NbSampEpoch;
          epochdef(iEpoch,3) = round(str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./(1000000./hdr.Fs)); % offset corresponds to timing
        end
      end
      
      if epochdef(end,2) ~= hdr.nSamples
        % check for NS 4.5.4 picosecond timing
        if (epochdef(end,2)/1000) == hdr.nSamples
          for iEpoch=1:size(epochdef,1)
            epochdef(iEpoch,1) = ((epochdef(iEpoch,1)-1)/1000)+1;
            epochdef(iEpoch,2) = epochdef(iEpoch,2)/1000;
            epochdef(iEpoch,3) = epochdef(iEpoch,3)/1000;
          end
          ft_warning('mff apparently generated by NetStation 4.5.4.  Adjusting time scale to microseconds from nanoseconds.');
        else
          ft_error('number of samples in all epochs do not add up to total number of samples')
        end
      end
      
      epochLengths = epochdef(:,2)-epochdef(:,1)+1;
      if ~any(diff(epochLengths))
        hdr.nSamples = epochLengths(1);
        hdr.nTrials  = length(epochLengths);
        
      else
        ft_warning('the data contains multiple epochs with variable length, possibly causing discontinuities in the data')
        % sanity check
        if epochdef(end,2) ~= hdr.nSamples
          % check for NS 4.5.4 picosecond timing
          if (epochdef(end,2)/1000) == hdr.nSamples
            for iEpoch=1:size(epochdef,1)
              epochdef(iEpoch,1)=((epochdef(iEpoch,1)-1)/1000)+1;
              epochdef(iEpoch,2)=epochdef(iEpoch,2)/1000;
              epochdef(iEpoch,3)=epochdef(iEpoch,3)/1000;
            end
            disp('mff apparently generated by NetStation 4.5.4.  Adjusting time scale to microseconds from nanoseconds.');
          else
            ft_error('number of samples in all epochs do not add up to total number of samples')
          end
        end
      end
      orig.epochdef = epochdef;
    end
    hdr.orig = orig;
    
  case 'egi_mff_v2'
    % ensure that the EGI_MFF_V2 toolbox is on the path
    ft_hastoolbox('egi_mff_v2', 1);
    
    %%%%%%%%%%%%%%%%%%%%%%
    %workaround for MATLAB bug resulting in global variables being cleared
    globalTemp=cell(0);
    globalList=whos('global');
    varList=whos;
    for i=1:length(globalList)
      eval(['global ' globalList(i).name ';']);
      eval(['globalTemp{end+1}=' globalList(i).name ';']);
    end
    %%%%%%%%%%%%%%%%%%%%%%
    
    % ensure that the JVM is running and the jar file is on the path
    mff_setup;
    
    %%%%%%%%%%%%%%%%%%%%%%
    %workaround for MATLAB bug resulting in global variables being cleared
    varNames={varList.name};
    for i=1:length(globalList)
      eval(['global ' globalList(i).name ';']);
      eval([globalList(i).name '=globalTemp{i};']);
      if ~any(strcmp(globalList(i).name,varNames)) %was global variable originally out of scope?
        eval(['clear ' globalList(i).name ';']); %clears link to global variable without affecting it
      end
    end
    clear globalTemp globalList varNames varList;
    %%%%%%%%%%%%%%%%%%%%%%
    
    if isunix && filename(1)~=filesep
      % add the full path to the dataset directory
      filename = fullfile(pwd, filename);
    elseif ispc && ~any(strcmp(filename(2),{':','\'}))
      % add the full path, including drive letter or slashes as needed.
      filename = fullfile(pwd, filename);
    end
    
    hdr = read_mff_header(filename);
    
  case {'egi_mff_v3' 'egi_mff'} % this is the default
    ft_hastoolbox('mffmatlabio', 1);
    hdr = mff_fileio_read_header(filename);
    
  case 'fcdc_buffer'
    % read from a networked buffer for realtime analysis
    [host, port] = filetype_check_uri(filename);
    
    if retry
      orig = [];
      while isempty(orig)
        try
          % try reading the header, catch the error and retry
          orig = buffer('get_hdr', [], host, port);
        catch
          ft_warning('could not read header from %s, retrying in 1 second', filename);
          pause(1);
        end
      end % while
    else
      % try reading the header only once, give error if it fails
      orig = buffer('get_hdr', [], host, port);
    end % if retry
    
    % construct the standard header elements
    hdr.Fs          = orig.fsample;
    hdr.nChans      = orig.nchans;
    hdr.nSamples    = orig.nsamples;
    hdr.nSamplesPre = 0;  % since continuous
    hdr.nTrials     = 1;  % since continuous
    hdr.orig        = []; % this will contain the chunks (if present)
    
    % add the contents of attached NEUROMAG_HEADER chunk after decoding to MATLAB structure
    if isfield(orig, 'neuromag_header')
      if isempty(cachechunk)
        % this only needs to be decoded once
        cachechunk = decode_fif(orig);
      end
      
      % convert to FieldTrip format header
      hdr.label       = cachechunk.ch_names(:);
      hdr.nChans      = cachechunk.nchan;
      hdr.Fs          = cachechunk.sfreq;
      
      % add a gradiometer structure for forward and inverse modelling
      try
        [grad, elec] = mne2grad(cachechunk, true, coilaccuracy); % the coordsys is 'dewar'
        if ~isempty(grad)
          hdr.grad = grad;
        end
        if ~isempty(elec)
          hdr.elec = elec;
        end
      catch
        disp(lasterr);
      end
      
      % store the original details
      hdr.orig = cachechunk;
    end
    
    % add the contents of attached CTF_RES4 chunk after decoding to MATLAB structure
    if isfield(orig, 'ctf_res4')
      if isempty(cachechunk)
        % this only needs to be decoded once
        cachechunk = decode_res4(orig.ctf_res4);
      end
      % copy the gradiometer details
      hdr.grad = cachechunk.grad;
      hdr.orig = cachechunk.orig;
      if isfield(orig, 'channel_names')
        % get the same selection of channels from the two chunks
        [selbuf, selres4] = match_str(orig.channel_names, cachechunk.label);
        if length(selres4)<length(orig.channel_names)
          ft_error('the res4 chunk did not contain all channels')
        end
        % copy some of the channel details
        hdr.label     = cachechunk.label(selres4);
        hdr.chantype  = cachechunk.chantype(selres4);
        hdr.chanunit  = cachechunk.chanunit(selres4);
        % add the channel names chunk as well
        hdr.orig.channel_names = orig.channel_names;
      end
      % add the raw chunk as well
      hdr.orig.ctf_res4 = orig.ctf_res4;
    end
    
    % add the contents of attached NIFTI_1 chunk after decoding to MATLAB structure
    if isfield(orig, 'nifti_1')
      hdr.nifti_1 = decode_nifti1(orig.nifti_1);
      % add the raw chunk as well
      hdr.orig.nifti_1 = orig.nifti_1;
    end
    
    % add the contents of attached SiemensAP chunk after decoding to MATLAB structure
    if isfield(orig, 'siemensap') && exist('sap2matlab')==3 % only run this if MEX file is present
      hdr.siemensap = sap2matlab(orig.siemensap);
      % add the raw chunk as well
      hdr.orig.siemensap = orig.siemensap;
    end
    
    if ~isfield(hdr, 'label')
      % prevent overwriting the labels that we might have gotten from a RES4 chunk
      if isfield(orig, 'channel_names')
        hdr.label = orig.channel_names;
      else
        hdr.label = cell(hdr.nChans,1);
        if hdr.nChans < 2000 % don't do this for fMRI etc.
          ft_warning('creating fake channel names');        % give this warning only once
          for i=1:hdr.nChans
            hdr.label{i} = sprintf('%d', i);
          end
        else
          ft_warning('skipping fake channel names');        % give this warning only once
          checkUniqueLabels = false;
        end
      end
    end
    
    if ~isfield(hdr, 'chantype')
      % prevent overwriting the chantypes that we might have gotten from a RES4 chunk
      hdr.chantype = cell(hdr.nChans,1);
      if hdr.nChans < 2000 % don't do this for fMRI etc.
        hdr.chantype = repmat({'unknown'}, 1, hdr.nChans);
      end
    end
    
    if ~isfield(hdr, 'chanunit')
      % prevent overwriting the chanunits that we might have gotten from a RES4 chunk
      hdr.chanunit = cell(hdr.nChans,1);
      if hdr.nChans < 2000 % don't do this for fMRI etc.
        hdr.chanunit = repmat({'unknown'}, 1, hdr.nChans);
      end
    end
    
    hdr.orig.bufsize = orig.bufsize;
    
    
  case 'fcdc_buffer_offline'
    [hdr, nameFlag] = read_buffer_offline_header(headerfile);
    switch nameFlag
      case 0
        % no labels generated (fMRI etc)
        checkUniqueLabels = false; % no need to check these
      case 1
        % has generated fake channels
        % give this warning only once
        ft_warning('creating fake channel names');
        checkUniqueLabels = false; % no need to check these
      case 2
        % got labels from chunk, check those
        checkUniqueLabels = true;
    end
    
  case 'fcdc_matbin'
    % this is multiplexed data in a *.bin file, accompanied by a MATLAB file containing the header
    load(headerfile, 'hdr');
    
  case 'fcdc_mysql'
    % check that the required low-level toolbox is available
    ft_hastoolbox('mysql', 1);
    % read from a MySQL server listening somewhere else on the network
    db_open(filename);
    if db_blob
      hdr = db_select_blob('fieldtrip.header', 'msg', 1);
    else
      hdr = db_select('fieldtrip.header', {'nChans', 'nSamples', 'nSamplesPre', 'Fs', 'label'}, 1);
      hdr.label = mxDeserialize(hdr.label);
    end
    
  case 'gtec_hdf5'
    % check that the required low-level toolbox is available
    ft_hastoolbox('gtec', 1);
    % there is only a precompiled *.p reader that reads the whole file at once
    orig = ghdf5read(filename);
    for i=1:numel(orig.RawData.AcquisitionTaskDescription.ChannelProperties.ChannelProperties)
      lab = orig.RawData.AcquisitionTaskDescription.ChannelProperties.ChannelProperties(i).ChannelName;
      typ = orig.RawData.AcquisitionTaskDescription.ChannelProperties.ChannelProperties(1).ChannelType;
      if isnumeric(lab)
        hdr.label{i} = num2str(lab);
      else
        hdr.label{i} = lab;
      end
      if ischar(typ)
        hdr.chantype{i} = lower(typ);
      else
        hdr.chantype{i} = 'unknown';
      end
    end
    hdr.Fs          = orig.RawData.AcquisitionTaskDescription.SamplingFrequency;
    hdr.nChans      = size(orig.RawData.Samples, 1);
    hdr.nSamples    = size(orig.RawData.Samples, 2);
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1; % assume continuous data, not epoched
    assert(orig.RawData.AcquisitionTaskDescription.NumberOfAcquiredChannels==hdr.nChans, 'inconsistent number of channels');
    % remember the complete data upon request
    if cache
      hdr.orig = orig;
    end
    
  case 'gtec_mat'
    % this is a simple MATLAB format, it contains a log and a names variable
    tmp = load(headerfile);
    log   = tmp.log;
    names = tmp.names;
    
    hdr.label = cellstr(names);
    hdr.nChans = size(log,1);
    hdr.nSamples = size(log,2);
    hdr.nSamplesPre = 0;
    hdr.nTrials = 1; % assume continuous data, not epoched
    
    % compute the sampling frequency from the time channel
    sel = strcmp(hdr.label, 'Time');
    time = log(sel,:);
    
    hdr.Fs = 1./(time(2)-time(1));
    
    % also remember the complete data upon request
    if cache
      hdr.orig.log = log;
      hdr.orig.names = names;
    end
    
  case 'gdf'
    % this requires the biosig toolbox
    ft_hastoolbox('BIOSIG', 1);
    % In the case that the gdf files are written by one of the FieldTrip
    % realtime applications, such as biosig2ft, the gdf recording can be
    % split over multiple 1GB files. The sequence of files is then
    %   filename.gdf   <- this is the one that should be specified as the filename/dataset
    %   filename_1.gdf
    %   filename_2.gdf
    %   ...
    
    [p, f, x] = fileparts(filename);
    if exist(sprintf('%s_%d%s', fullfile(p, f), 1, x), 'file')
      % there are multiple files, count the number of additional files (excluding the first one)
      count = 0;
      while exist(sprintf('%s_%d%s', fullfile(p, f), count+1, x), 'file')
        count = count+1;
      end
      hdr = read_biosig_header(filename);
      for i=1:count
        hdr(i+1) = read_biosig_header(sprintf('%s_%d%s', fullfile(p, f), i, x));
        % do some sanity checks
        if hdr(i+1).nChans~=hdr(1).nChans
          ft_error('multiple GDF files detected that should be appended, but the channel count is inconsistent');
        elseif hdr(i+1).Fs~=hdr(1).Fs
          ft_error('multiple GDF files detected that should be appended, but the sampling frequency is inconsistent');
        elseif ~isequal(hdr(i+1).label, hdr(1).label)
          ft_error('multiple GDF files detected that should be appended, but the channel names are inconsistent');
        end
      end % for count
      % combine all headers into one
      combinedhdr             = [];
      combinedhdr.Fs          = hdr(1).Fs;
      combinedhdr.nChans      = hdr(1).nChans;
      combinedhdr.nSamples    = sum([hdr.nSamples].*[hdr.nTrials]);
      combinedhdr.nSamplesPre = 0;
      combinedhdr.nTrials     = 1;
      combinedhdr.label       = hdr(1).label;
      combinedhdr.orig        = hdr; % include all individual file details
      hdr = combinedhdr;
      
    else
      % there is only a single file
      hdr = read_biosig_header(filename);
      % the GDF format is always continuous
      hdr.nSamples = hdr.nSamples * hdr.nTrials;
      hdr.nTrials = 1;
      hdr.nSamplesPre = 0;
    end % if single or multiple gdf files
    
  case {'homer_nirs'}
    % Homer files are MATLAB files in disguise
    % see https://www.nitrc.org/plugins/mwiki/index.php/homer2:Homer_Input_Files#NIRS_data_file_format
    nirs = load(filename, '-mat');
    % convert it to a raw data structure according to FT_DATATYPE_RAW
    data = homer2fieldtrip(nirs);
    % get the header information as structure
    hdr = ft_fetch_header(data);
    
  case {'itab_raw' 'itab_mhd'}
    % read the full header information frtom the binary header structure
    header_info = read_itab_mhd(headerfile);
    
    % these are the channels that are visible to FieldTrip
    chansel = 1:header_info.nchan;
    
    % convert the header information into a FieldTrip compatible format
    hdr.nChans      = length(chansel);
    hdr.label       = {header_info.ch(chansel).label};
    hdr.label       = hdr.label(:);  % should be column vector
    hdr.Fs          = header_info.smpfq;
    % it will always be continuous data
    hdr.nSamples    = header_info.ntpdata;
    hdr.nSamplesPre = 0; % it is a single continuous trial
    hdr.nTrials     = 1; % it is a single continuous trial
    % keep the original details AND the list of channels as used by FieldTrip
    hdr.orig         = header_info;
    hdr.orig.chansel = chansel;
    % add the gradiometer definition
    hdr.grad         = itab2grad(header_info);
    
  case 'jaga16'
    % this is hard-coded for the Jinga-Hi JAGA16 system with 16 channels
    packetsize = (4*2 + 6*2 + 16*43*2); % in bytes
    % read the first packet
    fid  = fopen_or_error(filename, 'r');
    buf  = fread(fid, packetsize/2, 'uint16');
    fclose(fid);
    
    if buf(1)==0
      % it does not have timestamps, i.e. it is the raw UDP stream
      packetsize = packetsize - 8; % in bytes
      packet     = jaga16_packet(buf(1:(packetsize/2)), false);
    else
      % each packet starts with a timestamp
      packet = jaga16_packet(buf, true);
    end
    
    % determine the number of packets from the file size
    info     = dir(filename);
    npackets = floor((info.bytes)/packetsize/2);
    
    hdr             = [];
    hdr.Fs          = packet.fsample;
    hdr.nChans      = packet.nchan;
    hdr.nSamples    = 43;
    hdr.nSamplesPre = 0;
    hdr.nTrials     = npackets;
    hdr.label       = cell(hdr.nChans,1);
    hdr.chantype    = cell(hdr.nChans,1);
    hdr.chanunit    = cell(hdr.nChans,1);
    for i=1:hdr.nChans
      hdr.label{i} = sprintf('%d', i);
      hdr.chantype{i} = 'eeg';
      hdr.chanunit{i} = 'uV';
    end
    
    % store some low-level details
    hdr.orig.offset     = 0;
    hdr.orig.packetsize = packetsize;
    hdr.orig.packet     = packet;
    hdr.orig.info       = info;
    
  case {'manscan_mbi', 'manscan_mb2'}
    orig       = in_fopen_manscan(filename);
    hdr.Fs     = orig.prop.sfreq;
    hdr.nChans = numel(orig.channelmat.Channel);
    hdr.nTrials  = 1;
    if isfield(orig, 'epochs') && ~isempty(orig.epochs)
      hdr.nSamples = 0;
      for i = 1:numel(orig.epochs)
        hdr.nSamples =  hdr.nSamples + diff(orig.epochs(i).samples) + 1;
      end
    else
      hdr.nSamples = diff(orig.prop.samples) + 1;
    end
    if orig.prop.times(1) < 0
      hdr.nSamplesPre  = round(orig.prop.times(1)/hdr.Fs);
    else
      hdr.nSamplesPre  = 0;
    end
    for i=1:hdr.nChans
      hdr.label{i,1}    = orig.channelmat.Channel(i).Name;
      hdr.chantype{i,1} = lower(orig.channelmat.Channel(i).Type);
      if isequal(hdr.chantype{i,1}, 'eeg')
        hdr.chanunit{i, 1} = 'uV';
      else
        hdr.chanunit{i, 1} = 'unknown';
      end
    end
    hdr.orig = orig;
    
  case 'matlab'
    % read the header structure from a MATLAB file
    % it should either contain a "hdr" structure, or a FieldTrip data structure according to FT_DATATYPE_RAW
    w = whos(matfile(filename));
    if any(strcmp({w.name}, 'hdr'))
      hdr = loadvar(filename, 'hdr');
    elseif any(strcmp({w.name}, 'data')) || length(w)==1
      data = loadvar(filename, 'data');
      hdr = ft_fetch_header(data);
    end

  case 'mayo_mef30'
    ft_hastoolbox('mayo_mef', 1); % make sure mayo_mef exists
    hdr = read_mayo_mef30(filename, password, sortchannel);
    
  case 'mayo_mef21'
    ft_hastoolbox('mayo_mef', 1); % make sure mayo_mef exists
    hdr = read_mayo_mef21(filename, password);
    
  case 'mega_neurone'
    % ensure that this external toolbox is on the path
    ft_hastoolbox('neurone', 1);
    if filename(end)~=filesep
      % it should end with a slash
      filename = [filename filesep];
    end
    % this is like the EEGLAB data structure
    EEG = readneurone(filename);
    
    hdr.Fs          = EEG.srate;
    hdr.nChans      = EEG.nbchan;
    hdr.nSamples    = EEG.pnts;
    hdr.nSamplesPre = -EEG.xmin*EEG.srate;
    hdr.nTrials     = EEG.trials;
    try
      hdr.label       = { EEG.chanlocs.labels }';
    catch
      ft_warning('creating default channel names');
      for i=1:hdr.nChans
        hdr.label{i} = sprintf('chan%03d', i);
      end
    end
    ind = 1;
    for i = 1:length( EEG.chanlocs )
      if isfield(EEG.chanlocs(i), 'X') && ~isempty(EEG.chanlocs(i).X)
        hdr.elec.label{ind, 1} = EEG.chanlocs(i).labels;
        % this channel has a position
        hdr.elec.elecpos(ind,1) = EEG.chanlocs(i).X;
        hdr.elec.elecpos(ind,2) = EEG.chanlocs(i).Y;
        hdr.elec.elecpos(ind,3) = EEG.chanlocs(i).Z;
        ind = ind+1;
      end
    end
    
    if cache
      % also remember the data and events
      hdr.orig = EEG;
    else
      % remember only the header details
      hdr.orig = removefields(EEG, {'data', 'event'});
    end
    
  case 'micromed_trc'
    orig = read_micromed_trc(filename);
    hdr             = [];
    hdr.Fs          = orig.Rate_Min; % FIXME is this correct?
    hdr.nChans      = orig.Num_Chan;
    hdr.nSamples    = orig.Num_Samples;
    hdr.nSamplesPre = 0; % continuous
    hdr.nTrials     = 1; % continuous
    hdr.label       = cell(1,hdr.nChans);
    % give this warning only once
    hdr.label  = {orig.elec.Name};
    hdr.chanunit = {orig.elec.Unit};
    hdr.subjectname = orig.name;
    %warning('using a modified read_micromed_trc() function');
    
    % this should be a column vector
    hdr.label = hdr.label(:);
    % remember the original header details
    hdr.orig = orig;
    
  case {'mpi_ds', 'mpi_dap'}
    hdr = read_mpi_ds(filename);
    
  case 'netmeg'
    ft_hastoolbox('netcdf', 1);
    
    % this will read all NetCDF data from the file and subsequently convert
    % each of the three elements into a more easy to parse MATLAB structure
    s = netcdf(filename);
    
    for i=1:numel(s.AttArray)
      fname = fixname(s.AttArray(i).Str);
      fval  = s.AttArray(i).Val;
      if ischar(fval)
        fval = fval(fval~=0); % remove the \0 characters
        fval = strtrim(fval); % remove insignificant whitespace
      end
      Att.(fname) = fval;
    end
    
    for i=1:numel(s.VarArray)
      fname = fixname(s.VarArray(i).Str);
      fval  = s.VarArray(i).Data;
      if ischar(fval)
        fval = fval(fval~=0); % remove the \0 characters
        fval = strtrim(fval); % remove insignificant whitespace
      end
      Var.(fname) = fval;
    end
    
    for i=1:numel(s.DimArray)
      fname = fixname(s.DimArray(i).Str);
      fval  = s.DimArray(i).Dim;
      if ischar(fval)
        fval = fval(fval~=0); % remove the \0 characters
        fval = strtrim(fval); % remove insignificant whitespace
      end
      Dim.(fname) = fval;
    end
    
    % convert the relevant fields into the default header structure
    hdr.Fs          = 1000/Var.samplinginterval;
    hdr.nChans      = length(Var.channelstatus);
    hdr.nSamples    = Var.numsamples;
    hdr.nSamplesPre = 0;
    hdr.nTrials     = size(Var.waveforms, 1);
    hdr.chanunit    = cellstr(reshape(Var.channelunits, hdr.nChans, 2));
    hdr.chantype    = cellstr(reshape(lower(Var.channeltypes), hdr.nChans, 3));
    
    ft_warning('creating fake channel names');
    hdr.label = cell(hdr.nChans, 1);
    for i=1:hdr.nChans
      hdr.label{i} = sprintf('%d', i);
    end
    
    % remember the original details of the file
    % note that this also includes the data
    % this is large, but can be reused elsewhere
    hdr.orig.Att = Att;
    hdr.orig.Var = Var;
    hdr.orig.Dim = Dim;
    
    % construct the gradiometer structure from the complete header information
    hdr.grad = netmeg2grad(hdr);
    
    
  case 'nervus_eeg'
    hdr = read_nervus_header(filename);
    checkUniqueLabels = false;
    
  case 'neuralynx_dma'
    hdr = read_neuralynx_dma(filename);
    
  case 'neuralynx_sdma'
    hdr = read_neuralynx_sdma(filename);
    
  case 'neuralynx_ncs'
    ncs = read_neuralynx_ncs(filename, 1, 0);
    [p, f, x]       = fileparts(filename);
    hdr.Fs          = ncs.hdr.SamplingFrequency;
    hdr.label       = {f};
    hdr.nChans      = 1;
    hdr.nTrials     = 1;
    hdr.nSamplesPre = 0;
    hdr.nSamples    = ncs.NRecords * 512;
    hdr.orig        = ncs.hdr;
    FirstTimeStamp  = ncs.hdr.FirstTimeStamp;  % this is the first timestamp of the first block
    LastTimeStamp   = ncs.hdr.LastTimeStamp;   % this is the first timestamp of the last block, i.e. not the timestamp of the last sample
    hdr.TimeStampPerSample = double(LastTimeStamp - FirstTimeStamp) ./ ((ncs.NRecords-1)*512);
    hdr.FirstTimeStamp     = FirstTimeStamp;
    
  case 'neuralynx_nse'
    nse = read_neuralynx_nse(filename, 1, 0);
    [p, f, x]       = fileparts(filename);
    hdr.Fs          = nse.hdr.SamplingFrequency;
    hdr.label       = {f};
    hdr.nChans      = 1;
    hdr.nTrials     = nse.NRecords;  % each record contains one waveform
    hdr.nSamples    = 32;            % there are 32 samples in each waveform
    hdr.nSamplesPre = 0;
    hdr.orig        = nse.hdr;
    % FIXME add hdr.FirstTimeStamp and hdr.TimeStampPerSample
    
  case {'neuralynx_ttl', 'neuralynx_tsl', 'neuralynx_tsh'}
    % these are hardcoded, they contain an 8-byte header and int32 values for a single channel
    % FIXME this should be done similar as neuralynx_bin, i.e. move the hdr into the function
    hdr             = [];
    hdr.Fs          = 32556;
    hdr.nChans      = 1;
    hdr.nSamples    = (filesize(filename)-8)/4;
    hdr.nSamplesPre = 1;
    hdr.nTrials     = 1;
    hdr.label       = {headerformat((end-3):end)};
    
  case 'neuralynx_bin'
    hdr = read_neuralynx_bin(filename);
    
  case 'neuralynx_ds'
    hdr = read_neuralynx_ds(filename);
    
  case 'neuralynx_cds'
    hdr = read_neuralynx_cds(filename);
    
  case 'nexstim_nxe'
    hdr = read_nexstim_nxe(filename);
    
  case 'neuromag_headpos'
    % neuromag headposition file created with maxfilter, the position varies over time
    orig = read_neuromag_headpos(filename);
    hdr.Fs          = 1/(orig.data(1,2)-orig.data(1,1));
    hdr.nChans      = size(orig.data,1);
    hdr.nSamples    = size(orig.data,2);
    hdr.nSamplesPre = 0;   % convert from ms to samples
    hdr.nTrials     = 1;
    hdr.label       = orig.colheaders;
    
  case 'neuromag_maxfilterlog'
    % neuromag log file created with maxfilter
    log = read_neuromag_maxfilterlog(filename);
    hdr = [];
    hdr.label = {'t' 'e' 'g' 'v' 'r' 'd'};
    for i=1:numel(log.hpi)
      for j=1:11
        hdr.label{end+1} = sprintf('hpi%d_%02d', i, j);
      end
    end
    hdr.nChans = length(hdr.label);
    hdr.nSamples = length(log.t);
    hdr.nSamplesPre = 0;
    hdr.nTrials = 1;
    hdr.Fs = 1 / median(diff(log.t));
    hdr.orig = log;
    
  case {'neuromag_fif' 'neuromag_mne'}
    % check that the required low-level toolbox is available
    ft_hastoolbox('mne', 1);
    
    info = fiff_read_meas_info(filename);
    
    % convert to FieldTrip format header
    hdr.label       = info.ch_names(:);
    hdr.nChans      = info.nchan;
    hdr.Fs          = info.sfreq;
    
    % add a gradiometer structure for forward and inverse modelling
    try
      [grad, elec] = mne2grad(info, strcmp(coordsys, 'dewar'), coilaccuracy);
      if ~isempty(grad)
        hdr.grad = grad;
      end
      if ~isempty(elec)
        hdr.elec = elec;
      end
    catch
      disp(lasterr);
    end
    
    iscontinuous  = 0;
    isepoched     = 0;
    isaverage     = 0;
    
    if isempty(fiff_find_evoked(filename)) % true if file contains no evoked responses
      try
        epochs = fiff_read_epochs(filename);
        isepoched = 1;
      catch
        % the "catch me" syntax is broken on MATLAB74, this fixes it
        me = lasterror;
        if strcmp(me.identifier, 'MNE:fiff_read_events')
          iscontinuous = 1;
        else
          rethrow(me)
        end
      end
      
    else
      isaverage = 1;
    end
    
    if iscontinuous
      try
        % we only use 1 input argument here to allow backward
        % compatibility up to MNE 2.6.x:
        raw = fiff_setup_read_raw(filename);
      catch
        % the "catch me" syntax is broken on MATLAB74, this fixes it
        me = lasterror;
        % there is an error - we try to use MNE 2.7.x (if present) to
        % determine if the cause is maxshielding:
        try
          allow_maxshield = true;
          raw = fiff_setup_read_raw(filename,allow_maxshield);
        catch
          % unknown problem, or MNE version 2.6.x or less:
          rethrow(me);
        end
        % no error message from fiff_setup_read_raw? Then maxshield
        % was applied, but maxfilter wasn't, so return this error:
        if istrue(checkmaxfilter)
          ft_error('Maxshield data should be corrected using Maxfilter prior to importing in FieldTrip.');
        else
          ft_warning('Maxshield data should be corrected using Maxfilter prior to importing in FieldTrip.');
        end
      end
      hdr.nSamples    = raw.last_samp - raw.first_samp + 1; % number of samples per trial
      hdr.nSamplesPre = 0;
      % otherwise conflicts will occur in read_data
      hdr.nTrials     = 1;
      info.raw        = raw; % keep all the details
      
    elseif isepoched
      hdr.nSamples    = length(epochs.times);
      hdr.nSamplesPre = sum(epochs.times < 0);
      hdr.nTrials     = size(epochs.data, 1);
      info.epochs     = epochs;  % this is used by read_data to get the actual data, i.e. to prevent re-reading
      
    elseif isaverage
      try
        evoked_data    = fiff_read_evoked_all(filename);
        vartriallength = any(diff([evoked_data.evoked.first])) || any(diff([evoked_data.evoked.last]));
        if vartriallength
          % there are trials averages with variable durations in the file
          ft_warning('EVOKED FILE with VARIABLE TRIAL LENGTH! - check data have been processed accurately');
          hdr.nSamples = 0;
          for i=1:length(evoked_data.evoked)
            hdr.nSamples = hdr.nSamples + size(evoked_data.evoked(i).epochs, 2);
          end
          % represent it as a continuous file with a single trial
          % all trial average details will be available through read_event
          hdr.nSamplesPre = 0;
          hdr.nTrials     = 1;
          info.evoked     = evoked_data.evoked; % this is used by read_data to get the actual data, i.e. to prevent re-reading
          info.info       = evoked_data.info;   % keep all the details
          info.vartriallength = 1;
        else
          % represent it as a file with multiple trials, each trial has the same length
          % all trial average details will be available through read_event
          hdr.nSamples    = evoked_data.evoked(1).last - evoked_data.evoked(1).first + 1;
          hdr.nSamplesPre = -evoked_data.evoked(1).first;   % represented as negative number in fif file
          hdr.nTrials     = length(evoked_data.evoked);
          info.evoked     = evoked_data.evoked;             % this is used by read_data to get the actual data, i.e. to prevent re-reading
          info.info       = evoked_data.info;               % keep all the details
          info.vartriallength = 0;
        end
      catch
        % this happens if fiff_read_evoked_all cannot find evoked
        % responses, in which case it errors due to not assigning the
        % output variable "data"
        ft_warning('%s does not contain data', filename);
        hdr.nSamples    = 0;
        hdr.nSamplesPre = 0;
        hdr.nTrials     = 0;
      end
    end
    
    % remember the original header details
    hdr.orig = info;
    
    % these are useful to know in ft_read_event and ft_read_data
    hdr.orig.isaverage    = isaverage;
    hdr.orig.iscontinuous = iscontinuous;
    hdr.orig.isepoched    = isepoched;
    
  case 'neuromag_mex'
    % check that the required low-level toolbox is available
    ft_hastoolbox('meg-pd', 1);
    rawdata('any',filename);
    rawdata('goto', 0);
    megmodel('head',[0 0 0],filename);
    % get the available information from the fif file
    [orig.rawdata.range,orig.rawdata.calib]           = rawdata('range');
    [orig.rawdata.sf]                                 = rawdata('sf');
    [orig.rawdata.samples]                            = rawdata('samples');
    [orig.chaninfo.N,orig.chaninfo.S,orig.chaninfo.T] = chaninfo;           % Numbers, names & places
    [orig.chaninfo.TY,orig.chaninfo.NA]               = chaninfo('type');   % Coil type
    [orig.chaninfo.NO]                                = chaninfo('noise');  % Default noise level
    [orig.channames.NA,orig.channames.KI,orig.channames.NU] = channames(filename); % names, kind, logical numbers
    % read a single trial to determine the data size
    [buf, status] = rawdata('next');
    rawdata('close');
    
    % This is to solve a problem reported by Doug Davidson: The problem
    % is that rawdata('samples') is not returning the number of samples
    % correctly. It appears that the example script rawchannels in meg-pd
    % might work, however, so I want to use rawchannels to read in one
    % channel of data in order to get the number of samples in the file:
    if orig.rawdata.samples<0
      tmpchannel = 1;
      tmpvar = rawchannels(filename,tmpchannel);
      [orig.rawdata.samples] = size(tmpvar,2);
      clear tmpvar tmpchannel;
    end
    
    % convert to FieldTrip format header
    hdr.label       = orig.channames.NA;
    hdr.Fs          = orig.rawdata.sf;
    hdr.nSamplesPre = 0; % I don't know how to get this out of the file
    hdr.nChans      = size(buf,1);
    hdr.nSamples    = size(buf,2); % number of samples per trial
    hdr.nTrials     = orig.rawdata.samples ./ hdr.nSamples;
    % add a gradiometer structure for forward and inverse modelling
    hdr.grad = fif2grad(filename);
    % remember the original header details
    hdr.orig = orig;
    
  case 'neuroomega_mat'
    % These are MATLAB *.mat files created by the software 'Map File
    % Converter' from the proprietary .mpx files recorded by NeuroOmega
    
    %defining default if chantype is missing
    if isempty(chantype) || strcmpi(chantype{1},'chaninfo') || strcmpi(chantype{1},'info')
      chantype={'micro'};
    end
    
    chantype_dict={'micro','macro',     'analog', 'micro_lfp','macro_lfp','micro_hp','add_analog','emg', 'eeg';...
      'CRAW', 'CMacro_RAW','CANALOG', 'CLFP',     'CMacro_LFP',   'CSPK' ,'CADD_ANALOG','CEMG', 'CEEG'};
    neuroomega_param={'_KHz','_KHz_Orig','_Gain','_BitResolution','_TimeBegin','_TimeEnd'};
    
    %identifying channels to be loaded
    orig = matfile(filename);
    fields_orig=who(orig);
    
    %is_param=endsWith(fields_orig,neuroomega_param); %Matlab 2017a
    is_param=zeros(length(fields_orig),1);
    for i=1:length(neuroomega_param)
      is_param = is_param | ~cellfun('isempty',regexp(fields_orig,strcat(neuroomega_param(i),'$'),'start'));
    end
    
    %creating channel info table
    channels_all=fields_orig(contains(fields_orig,chantype_dict(2,:)) & ~is_param);
    %Matching channels to chantypes
    M=cell2mat(cellfun(@(x) contains(channels_all,x),chantype_dict(2,:),'UniformOutput',false));
    chantype_ix = sum( cumprod(M == 0, 2), 2) + 1;
    Fs=cellfun(@(x) orig.([x '_KHz'])*1000,channels_all);
    chaninfo=table(channels_all,chantype_dict(1,chantype_ix)',Fs,'VariableNames',{'channel' 'chantype' 'Fs'});
    
    channels={}; channelstype={};
    for c = 1:length(chantype)
      chantype_dict_sel=strcmpi(chantype_dict(1,:),chantype{c});
      if sum(chantype_dict_sel)>0
        chanbasename=chantype_dict{2, chantype_dict_sel};
        sel_chan=fields_orig(strncmpi(fields_orig,chanbasename,length(chanbasename)) & ~is_param);
        if isempty(sel_chan)
          ft_warning(strjoin({'chantype ',chantype{c},' selected but no ',chanbasename,' found'}))
        else
          channels=[channels;sel_chan];
          channelstype=[channelstype;repmat(chantype(c),  size(sel_chan))];
        end
      else
        ft_warning(strjoin({'unknown chantype ',chantype{c}}));
      end
    end
    
    if ~isempty(channels)
      chan_t=table;
      for i=1:length(channels)
        ch=channels{i};
        ch_whos=whos(orig,ch);
        chan_t=[chan_t;{ch,orig.([ch,'_KHz'])*1000, orig.([ch,'_TimeBegin']), ch_whos.size(2)}];
      end
      chan_t.Properties.VariableNames={'channel' 'Fs' 'T0' 'nSamples'};
      
      Fs=unique(chan_t.Fs);
      T0=unique(chan_t.T0);
      nSamples=unique(chan_t.nSamples);
      
      if length(Fs)>1
        chan_t %; printing table for user
        ft_error('inconsistent channels with different sampling rates for %s',filename);
      end
      if length(T0)>1
        chan_t %; printing table for user
        ft_warning('inconsistent channels with different initial times for %s. Selecting minimum time',filename);
        T0 = min(T0);
      end
      if length(nSamples)>1
        chan_t %; printing table for user
        ft_warning('inconsistent number of samples across channels for %s. Selecting minimun nSample',filename)
        nSamples=min(nSamples);
      end
    else %If no channel selected
      Fs=nan; nSamples=nan; channelstype=chaninfo.chantype;
    end
    
    % building header
    hdr.Fs          = Fs;
    hdr.nChans      = length(channels);
    hdr.nSamples    = nSamples;
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1;
    hdr.label       = deblank(channels);
    hdr.chantype    = channelstype;
    hdr.chanunit    = repmat({'uV'},  size(hdr.label));
    % store the complete information in hdr.orig
    % ft_read_data and ft_read_event will get it from there
    hdr.orig        = [];
    hdr.orig.orig   = orig;
    hdr.orig.chaninfo = chaninfo;
    hdr.orig.fields = fields_orig;
    
  case 'neuroprax_eeg'
    orig = np_readfileinfo(filename);
    
    hdr.Fs          = orig.fa;
    hdr.nChans      = orig.K;
    hdr.nSamples    = orig.N;
    hdr.nSamplesPre = 0; % continuous
    hdr.nTrials     = 1; % continuous
    hdr.label       = orig.channels(:);
    hdr.unit        = orig.units(:);
    
    % remember the original header details
    hdr.orig = orig;
    
  case 'neuroscope_bin'
    [p,f,e]    = fileparts(filename);
    headerfile = fullfile(p,[f,'.xml']);
    hdr        = ft_read_header(headerfile, 'headerformat', 'neuroscope_xml');
    
  case 'neuroscope_ds'
    listing    = dir(filename);
    filenames  = {listing.name}';
    headerfile = filenames{~cellfun('isempty',strfind(filenames,'.xml'))};
    hdr        = ft_read_header(headerfile, 'headerformat', 'neuroscope_xml');
    
  case 'neuroscope_xml'
    ft_hastoolbox('neuroscope', 1);
    ft_hastoolbox('gifti', 1);
    
    % this pertains to generic header file, and the other neuroscope
    % formats will recurse into this one
    [p,f,e]    = fileparts(filename);
    if isempty(p), p = pwd; end
    listing    = dir(p);
    filenames  = {listing.name}';
    
    lfpfile_idx = find(~cellfun('isempty',strfind(filenames,'.eeg')));
    rawfile_idx = find(~cellfun('isempty',strfind(filenames,'.dat')));
    
    if ~isempty(lfpfile_idx)
      % FIXME this assumes only 1 such file, or at least it only takes the
      % first one.
      lfpfile = filenames{lfpfile_idx(1)};
    end
    if ~isempty(rawfile_idx)
      rawfile = filenames{rawfile_idx(1)};
    end
    params     = LoadParameters(filename);
    
    hdr         = [];
    hdr.nChans  = params.nChannels;
    hdr.nTrials = 1; % is it always continuous? FIXME
    hdr.nSamplesPre = 0;
    
    if ~isempty(lfpfile)
      % use the sampling of the lfp-file to be leading
      hdr.Fs       = params.rates.lfp;
      hdr.nSamples = listing(strcmp(filenames,lfpfile)).bytes./(hdr.nChans*params.nBits/8);
      hdr.TimeStampPerSample = params.rates.wideband./params.rates.lfp;
    else
      % use the sampling of the raw-file to be leading
      hdr.Fs     = params.rates.wideband;
      hdr.nSamples = listing(strcmp(filenames,rawfile)).bytes./(hdr.nChans*params.nBits/8);
      hdr.TimeStampPerSample = 1;
    end
    hdr.orig = params;
    
    hdr.label = cell(hdr.nChans,1);
    for k = 1:hdr.nChans
      hdr.label{k} = ['chan',num2str(k,'%0.3d')];
    end
    
  case 'neurosim_evolution'
    hdr = read_neurosim_evolution(filename);
    
  case {'neurosim_ds' 'neurosim_signals'}
    hdr = read_neurosim_signals(filename);
    
  case 'neurosim_spikes'
    headerOnly = true;
    hdr = read_neurosim_spikes(filename, headerOnly);
    
  case 'nihonkohden_m00'
    % this is an ASCII file format which is rather inefficient to read
    if cache
      % read it once and store the data along with the header
      [hdr, dat] = read_nihonkohden_m00(filename);
      hdr.orig.dat = dat;
    else
      % read only the header
      hdr = read_nihonkohden_m00(filename);
    end
    
  case 'nihonkohden_eeg'
    ft_hastoolbox('brainstorm', 1);
    hdr = read_brainstorm_header(filename);
    
  case 'nimh_cortex'
    cortex = read_nimh_cortex(filename, 'epp', 'no', 'eog', 'no');
    % look at the first trial to determine whether it contains data in the EPP and EOG channels
    trial1  = read_nimh_cortex(filename, 'epp', 'yes', 'eog', 'yes', 'begtrial', 1, 'endtrial', 1);
    hasepp = ~isempty(trial1.epp);
    haseog = ~isempty(trial1.eog);
    if hasepp
      ft_warning('EPP channels are not yet supported');
    end
    % at the moment only the EOG channels are supported here
    if haseog
      hdr.label       = {'EOGx' 'EOGy'};
      hdr.nChans      = 2;
    else
      hdr.label       = {};
      hdr.nChans      = 0;
    end
    hdr.nTrials     = length(cortex);
    hdr.nSamples    = inf;
    hdr.nSamplesPre = 0;
    hdr.orig.trial = cortex;
    hdr.orig.hasepp = hasepp;
    hdr.orig.haseog = haseog;
    
  case 'ns_avg'
    orig = read_ns_hdr(filename);
    % do some reformatting/renaming of the header items
    hdr.Fs          = orig.rate;
    hdr.nSamples    = orig.npnt;
    hdr.nSamplesPre = round(-orig.rate*orig.xmin/1000);
    hdr.nChans      = orig.nchan;
    hdr.label       = orig.label(:);
    hdr.nTrials     = 1; % the number of trials in this datafile is only one, i.e. the average
    % remember the original header details
    hdr.orig = orig;
    
  case {'ns_cnt' 'ns_cnt16', 'ns_cnt32'}
    ft_hastoolbox('eeglab', 1);
    if strcmp(headerformat, 'ns_cnt')
      orig = loadcnt(filename); % let loadcnt figure it out
    elseif strcmp(headerformat, 'ns_cnt16')
      orig = loadcnt(filename, 'dataformat', 'int16');
    elseif strcmp(headerformat, 'ns_cnt32')
      orig = loadcnt(filename, 'dataformat', 'int32');
    end
    
    % do some reformatting/renaming of the header items
    hdr.Fs          = orig.header.rate;
    hdr.nChans      = orig.header.nchannels;
    hdr.nSamples    = orig.ldnsamples;
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1;
    for i=1:hdr.nChans
      hdr.label{i} = deblank(orig.electloc(i).lab);
    end
    % remember the original header details
    hdr.orig = orig;
    
  case 'ns_eeg'
    orig = read_ns_hdr(filename);
    % do some reformatting/renaming of the header items
    hdr.label       = orig.label;
    hdr.Fs          = orig.rate;
    hdr.nSamples    = orig.npnt;
    hdr.nSamplesPre = round(-orig.rate*orig.xmin/1000);
    hdr.nChans      = orig.nchan;
    hdr.nTrials     = orig.nsweeps;
    % remember the original header details
    hdr.orig = orig;
    
  case 'nmc_archive_k'
    hdr = read_nmc_archive_k_hdr(filename);
    
  case 'neuroshare' % NOTE: still under development
    % check that the required neuroshare toolbox is available
    ft_hastoolbox('neuroshare', 1);
    tmp = read_neuroshare(filename);
    hdr.Fs          = tmp.hdr.analoginfo(end).SampleRate; % take the sampling freq from the last analog channel (assuming this is the same for all chans)
    hdr.nChans      = length(tmp.list.analog(tmp.analog.contcount~=0)); % get the analog channels, only the ones that are not empty
    hdr.nSamples    = max([tmp.hdr.entityinfo(tmp.list.analog).ItemCount]); % take the number of samples from the longest channel
    hdr.nSamplesPre = 0; % continuous data
    hdr.nTrials     = 1; % continuous data
    hdr.label       = {tmp.hdr.entityinfo(tmp.list.analog(tmp.analog.contcount~=0)).EntityLabel}; %%% contains non-unique chans?
    hdr.orig        = tmp; % remember the original header
    
  case 'nwb'
    ft_hastoolbox('MatNWB', 1);	% when I run this locally outside of ft_read_header it does not work for me
    try
      c = load('namespaces/core.mat');
      nwb_version = c.version;
      nwb_fileversion = util.getSchemaVersion(filename);
      if ~strcmp(nwb_version, nwb_fileversion)
        warning(['Installed NWB:N schema version (' nwb_version ') does not match the file''s schema (' nwb_fileversion{1} '). This might result in an error. If so, try to install the matching schema from here: https://github.com/NeurodataWithoutBorders/nwb-schema/releases'])
      end
    catch
      warning('Something might not be alright with your MatNWB path. Will try anyways.')
    end
    tmp = nwbRead(filename); % is lazy, so should not be too costly
    es_key = tmp.searchFor('ElectricalSeries').keys; % find lfp data, which should be an ElectricalSeries object
    es_key = es_key(~contains(es_key, 'acquisition'));
    if isempty(es_key)
      error('Dataset does not contain an LFP signal (i.e., no object of the class ''ElectricalSeries''.')
    elseif numel(es_key) > 1 % && isempty(additional_user_input) % TODO: Try to sort this out with the user's help
      % Temporary fix: SpikeEventSeries is a daughter of ElectrialSeries but should not be found here (searchFor update on its way)
      es_key = es_key(contains(es_key,'lfp','IgnoreCase',true));
    end
    if numel(es_key) > 1 % in case we weren't able to sort out a single
      error('More than one ElectricalSeries present in data. Please specify which signal to use.')
    else
      eseries = io.resolvePath(tmp, es_key{1});
    end
    if isa(eseries.data, 'types.untyped.DataStub')
      hdr.nSamples = eseries.data.dims(2);
    elseif isa(eseries.data, 'types.untyped.DataPipe')
      hdr.nSamples = eseries.data.internal.maxSize(2);
    else
      warning('Cannot determine number of samples in the data.')
      hdr.nSamples = [];
    end
    hdr.Fs          = eseries.starting_time_rate;
    hdr.nSamplesPre = 0; % for now: hardcoded continuous data
    hdr.nTrials     = 1; % for now: hardcoded continuous data
    hdr.label       = {};
    tmp_ch          = io.resolvePath(tmp, eseries.electrodes.table.path).id.data.load; % electrode names
    for iCh=1:numel(tmp_ch) % TODO: does that work if nwb ids are strings?
      if isnumeric(tmp_ch(iCh))
        hdr.label(iCh,1) = {num2str(tmp_ch(iCh))};
      else
        hdr.label(iCh,1) = tmp_ch(iCh);
      end
    end
    hdr.nChans      = numel(hdr.label);
    [hdr.chanunit{1:hdr.nChans,1}] = deal(eseries.data_unit);
    hdr.chanunit    = strrep(hdr.chanunit, 'volt', 'V');
    hdr.chanunit    = strrep(hdr.chanunit, 'micro', 'u');
    % TODO: hdr.FirstTimeStamp
    % TODO: hdr.TimeStampPerSample
    
    % carry over some metadata
    hdr.orig        = [];
    fn = {'general_experimenter', ...
      'general_institution', ...
      'general_keywords', ...
      'general_lab', ...
      'general_notes', ...
      'general_related_publications', ...
      'general_session_id', ...
      'identifier', ...
      'session_description', ...
      'nwb_version', ...
      'help'};
    for iFn = 1:numel(fn)
      if isprop(tmp, fn{iFn}) && ~isempty(tmp.(fn{iFn}))
        hdr.orig.(fn{iFn}) = tmp.(fn{iFn});
      end
    end
    
  case 'artinis_oxy3'
    ft_hastoolbox('artinis', 1);
    hdr = read_artinis_oxy3(filename);
    
  case 'artinis_oxy4'
    ft_hastoolbox('artinis', 1);
    hdr = read_artinis_oxy4(filename);
    
  case 'artinis_oxyproj'
    ft_hastoolbox('artinis', 1);
    hdr = read_oxyproj_header(filename);
    
  case 'plexon_ds'
    hdr = read_plexon_ds(filename);
    
  case 'plexon_ddt'
    orig = read_plexon_ddt(filename);
    hdr.nChans      = orig.NChannels;
    hdr.Fs          = orig.Freq;
    hdr.nSamples    = orig.NSamples;
    hdr.nSamplesPre = 0;      % continuous
    hdr.nTrials     = 1;      % continuous
    hdr.label       = cell(1,hdr.nChans);
    % give this warning only once
    ft_warning('creating fake channel names');
    for i=1:hdr.nChans
      hdr.label{i} = sprintf('%d', i);
    end
    % also remember the original header
    hdr.orig        = orig;
    
  case {'read_nex_data'} % this is an alternative reader for nex files
    orig = read_nex_header(filename);
    % assign the obligatory items to the output FCDC header
    numsmp = cell2mat({orig.varheader.numsmp});
    adindx = find(cell2mat({orig.varheader.typ})==5);
    if isempty(adindx)
      ft_error('file does not contain continuous channels');
    end
    hdr.nChans      = length(orig.varheader);
    hdr.Fs          = orig.varheader(adindx(1)).wfrequency;     % take the sampling frequency from the first A/D channel
    hdr.nSamples    = max(numsmp(adindx));                      % take the number of samples from the longest A/D channel
    hdr.nTrials     = 1;                                        % it can always be interpreted as continuous data
    hdr.nSamplesPre = 0;                                        % and therefore it is not trial based
    for i=1:hdr.nChans
      hdr.label{i} = deblank(char(orig.varheader(i).nam));
    end
    hdr.label = hdr.label(:);
    % also remember the original header details
    hdr.orig = orig;
    
  case {'plexon_nex' 'read_plexon_nex'} % this is the default reader for nex files
    orig = read_plexon_nex(filename);
    numsmp = cell2mat({orig.VarHeader.NPointsWave});
    adindx = find(cell2mat({orig.VarHeader.Type})==5);
    if isempty(adindx)
      ft_error('file does not contain continuous channels');
    end
    hdr.nChans      = length(orig.VarHeader);
    hdr.Fs          = orig.VarHeader(adindx(1)).WFrequency;     % take the sampling frequency from the first A/D channel
    hdr.nSamples    = max(numsmp(adindx));                      % take the number of samples from the longest A/D channel
    hdr.nTrials     = 1;                                        % it can always be interpreted as continuous data
    hdr.nSamplesPre = 0;                                        % and therefore it is not trial based
    for i=1:hdr.nChans
      hdr.label{i} = deblank(char(orig.VarHeader(i).Name));
    end
    hdr.label = hdr.label(:);
    hdr.FirstTimeStamp     = orig.FileHeader.Beg;
    hdr.TimeStampPerSample = orig.FileHeader.Frequency ./ hdr.Fs;
    % also remember the original header details
    hdr.orig = orig;
    
  case 'plexon_nex5' % this is the default reader for nex5 files
    orig = read_nex5(filename);
    numsmp = cell2mat({orig.VarHeader.NumberOfDataPoints});
    adindx = find(cell2mat({orig.VarHeader.Type})==5);
    if isempty(adindx)
      ft_error('file does not contain continuous channels');
    end
    % check that all continuous channels have the same sampling rate
    samplingRates = cell2mat({orig.VarHeader.WFrequency});
    contSamplingRates = samplingRates(adindx);
    if any(contSamplingRates~=contSamplingRates(1))
      ft_error('different sampling rates in continuous data not supported');
    end
    hdr.nChans      = length(orig.VarHeader);
    hdr.Fs          = orig.VarHeader(adindx(1)).WFrequency;    % take the sampling frequency from the first A/D channel
    hdr.TimeStampPerSample = orig.FileHeader.Frequency ./ hdr.Fs;
    % for hdr.nSamples, we need to calculate the last timestamp for every continuous channel
    maxTimestamp = 0;
    for i = 1:length(adindx)
      [nex, chanhdr] = read_nex5(filename, 'header', orig, 'channel', adindx(i), 'tsonly', 1);
      numPointsInLastFragment = numsmp(adindx(i)) - nex.indx(end) - 1;
      maxTimestamp = max(maxTimestamp, nex.ts(end)+hdr.TimeStampPerSample*(numPointsInLastFragment-1));
    end
    
    hdr.nSamples    = maxTimestamp/hdr.TimeStampPerSample;
    hdr.nTrials     = 1;                                        % it can always be interpreted as continuous data
    hdr.nSamplesPre = 0;                                        % and therefore it is not trial based
    for i=1:hdr.nChans
      hdr.label{i} = deblank(char(orig.VarHeader(i).Name));
    end
    hdr.label = hdr.label(:);
    hdr.FirstTimeStamp     = orig.FileHeader.Beg;
    hdr.TimeStampPerSample = orig.FileHeader.Frequency ./ hdr.Fs;
    % also remember the original header details
    hdr.orig = orig;
    
  case 'plexon_plx'
    orig = read_plexon_plx(filename);
    if orig.NumSlowChannels==0
      ft_error('file does not contain continuous channels');
    end
    fsample = [orig.SlowChannelHeader.ADFreq];
    if any(fsample~=fsample(1))
      ft_error('different sampling rates in continuous data not supported');
    end
    for i=1:length(orig.SlowChannelHeader)
      label{i} = deblank(orig.SlowChannelHeader(i).Name);
    end
    % continuous channels don't always contain data, remove the empty ones
    sel  = [orig.DataBlockHeader.Type]==5;  % continuous
    chan = [orig.DataBlockHeader.Channel];
    for i=1:length(label)
      chansel(i) = any(chan(sel)==orig.SlowChannelHeader(i).Channel);
    end
    chansel = find(chansel); % this is required for timestamp selection
    label = label(chansel);
    % only the continuous channels are returned as visible
    hdr.nChans      = length(label);
    hdr.Fs          = fsample(1);
    hdr.label       = label;
    % also remember the original header
    hdr.orig        = orig;
    
    % select the first continuous channel that has data
    sel = ([orig.DataBlockHeader.Type]==5 & [orig.DataBlockHeader.Channel]==orig.SlowChannelHeader(chansel(1)).Channel);
    % get the timestamps that correspond with the continuous data
    tsl = [orig.DataBlockHeader(sel).TimeStamp]';
    tsh = [orig.DataBlockHeader(sel).UpperByteOf5ByteTimestamp]';
    ts  = timestamp_plexon(tsl, tsh);  % use helper function, this returns an uint64 array
    
    % determine the number of samples in the continuous channels
    num = [orig.DataBlockHeader(sel).NumberOfWordsInWaveform];
    hdr.nSamples    = sum(num);
    hdr.nSamplesPre = 0;      % continuous
    hdr.nTrials     = 1;      % continuous
    
    % the timestamps indicate the beginning of each block, hence the timestamp of the last block corresponds with the end of the previous block
    hdr.TimeStampPerSample = double(ts(end)-ts(1))/sum(num(1:(end-1)));
    hdr.FirstTimeStamp     = ts(1);                                                %  the timestamp of the first continuous sample
    
    % also make the spike channels visible
    for i=1:length(orig.ChannelHeader)
      hdr.label{end+1} = deblank(orig.ChannelHeader(i).Name);
    end
    hdr.label = hdr.label(:);
    hdr.nChans = length(hdr.label);
    
  case {'ricoh_ave', 'ricoh_con', 'ricoh_mrk'}
    % header can be read with Ricoh MEG Reader
    hdr = read_ricoh_header(filename);
    % add a gradiometer structure for forward and inverse modelling
    hdr.grad = ricoh2grad(hdr);
    hdr.chantype = ft_chantype(hdr.label);
    unk = find(strcmp('unknown', hdr.chantype));
    %  Warning message:
    if ~isempty(unk)
      label_unk = hdr.label(unk);
      no_unk = num2cell(unk);
      C = [label_unk(:), no_unk(:)] .';
      ft_warning(['Unknown channel types: (label, no) =' repmat('( %s, %d ) ', 1, length(unk) ) '\n'], C{:})
    end
    
  case 'smi_txt'
    smi = read_smi_txt(filename);
    hdr.nChans              = size(smi.dat,1);
    hdr.nSamples            = size(smi.dat,2);
    hdr.nSamplesPre         = 0;
    hdr.nTrials             = 1;
    
    hdr.label               = smi.label;
    hdr.Fs                  = smi.Fs;
    hdr.FirstTimeStamp      = smi.timestamp(1);
    hdr.TimeStampPerSample  = mean(diff(smi.timestamp)); % these timestamps are in microseconds
    
    if cache
      % remember all header and data details upon request
      hdr.orig = smi;
    else
      % remember only the original header details
      hdr.orig.header = smi.header;
    end
    
    % add channel units when possible.
    for i=1:hdr.nChans
      chanunit = regexp(hdr.label{i,1},'(?<=\[).+?(?=\])','match');
      if ~isempty(chanunit)
        hdr.chanunit{i,1} = chanunit{1};
        hdr.chantype{i,1} = 'eyetracker';
      else
        hdr.chanunit{i,1} = 'unknown';
        hdr.chantype{i,1} = 'unknown';
      end
    end
    
  case 'tmsi_poly5'
    orig = read_tmsi_poly5(filename);
    % the header contains all channels twice (for the low and high data word)
    % it seems that the file format was designed for 16 bit, and only later extended to 32 bit
    hdr             = [];
    hdr.nChans      = orig.header.NumberOfSignals/2;
    hdr.Fs          = orig.header.FS;
    hdr.nSamples    = orig.header.NumberSampleBlocks * orig.header.SamplePeriodsPerBlock;
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1; % continuous
    for i=2:2:orig.header.NumberOfSignals
      % remove the '(Lo) ' and the '(Hi) ' section
      hdr.label{i/2} = strtrim(orig.description(i).SignalName(6:end)');
    end
    % determine the EEG channels
    iseeg = true(size(hdr.label));
    iseeg = iseeg & cellfun(@isempty, regexp(hdr.label, 'BIP.*'));
    iseeg = iseeg & cellfun(@isempty, regexp(hdr.label, 'AUX.*'));
    iseeg = iseeg & cellfun(@isempty, regexp(hdr.label, 'Digi.*'));
    iseeg = iseeg & cellfun(@isempty, regexp(hdr.label, 'Saw.*'));
    iseeg = iseeg & cellfun(@isempty, regexp(hdr.label, 'Bit.*'));
    istrg = ~cellfun(@isempty, regexp(hdr.label, 'Digi.*'));
    hdr.chanunit = cell(size(hdr.label));
    hdr.chantype = cell(size(hdr.label));
    hdr.chanunit(:) = {'unknown'};
    hdr.chantype(:) = {'unknown'};
    hdr.chanunit(iseeg) = {'uV'};
    hdr.chantype(iseeg) = {'eeg'};
    hdr.chantype(istrg) = {'trigger'};
    % remember the original header details
    hdr.orig = orig;
    
  case 'tobii_tsv'
    tsv = read_tobii_tsv(filename);
    % keyboard
    % remember the original header details
    hdr.orig = tsv;
    
  case {'tdt_tsq', 'tdt_tev'}
    % FIXME the code below is not yet functional, it requires more input from the ESI in Frankfurt
    %     tsq = read_tdt_tsq(headerfile);
    %     k = 0;
    %     chan = unique([tsq.channel]);
    %     % loop over the physical channels
    %     for i=1:length(chan)
    %       chansel = [tsq.channel]==chan(i);
    %       code = unique({tsq(chansel).code});
    %       % loop over the logical channels
    %       for j=1:length(code)
    %         codesel = false(size(tsq));
    %         for k=1:numel(codesel)
    %           codesel(k) = isequal(tsq(k).code, code{j});
    %         end
    %         % find the first instance of this logical channel
    %         this = find(chansel(:) & codesel(:), 1);
    %         % add it to the list of channels
    %         k = k + 1;
    %         frequency(k) = tsq(this).frequency;
    %         label{k}     = [char(typecast(tsq(this).code, 'uint8')) num2str(tsq(this).channel)];
    %         tsqorig(k)   = tsq(this);
    %       end
    %     end
    ft_error('not yet implemented');
    
  case {'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw', 'yokogawa_mrk'}
    % header can be read with two toolboxes: Yokogawa MEG Reader and Yokogawa MEG160 (old inofficial toolbox)
    % newest toolbox takes precedence.
    if ft_hastoolbox('yokogawa_meg_reader', 3); % stay silent if it cannot be added
      hdr = read_yokogawa_header_new(filename);
      % add a gradiometer structure for forward and inverse modelling
      hdr.grad = yokogawa2grad_new(hdr);
      hdr.chantype = ft_chantype(hdr.label);
      unk = find(strcmp('unknown', hdr.chantype));
      %  Warning message:
      if ~isempty(unk)
        label_unk = hdr.label(unk);
        no_unk = num2cell(unk);
        C = [label_unk(:), no_unk(:)] .';
        ft_warning(['Unknown channel types: (label, no) =' repmat('( %s, %d ) ', 1, length(unk) ) '\n'], C{:})
      end
    else
      ft_hastoolbox('yokogawa', 1); % try it with the old version of the toolbox
      hdr = read_yokogawa_header(filename);
      % add a gradiometer structure for forward and inverse modelling
      hdr.grad = yokogawa2grad(hdr);
    end
    
    
  case {'audio_wav', 'audio_ogg', 'audio_flac', 'audio_au', 'audio_aiff', 'audio_aif', 'audio_aifc', 'audio_mp3', 'audio_m4a', 'audio_mp4'}
    % prior to MATLAB R2015b this used to be done with "wavread"
    % but the audioinfo/audioread function are at least available from 2012b up
    info = audioinfo(filename);
    hdr.Fs          = info.SampleRate;
    hdr.nChans      = info.NumChannels;
    hdr.nSamples    = info.TotalSamples;
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1;
    [p, f, x] = fileparts(filename);
    if hdr.nChans>1
      for i=1:hdr.nChans
        % use the file name and channel number
        hdr.label{i,1} = sprintf('%s channel %d', f, i);
        hdr.chantype{i,1} = 'audio';
      end
    else
      hdr.label{1,1} = f;
      hdr.chantype{1,1} = 'audio';
    end
    % remember the details
    hdr.orig = info;
    
  case 'videomeg_aud'
    hdr = read_videomeg_aud(filename);
    
  case 'videomeg_vid'
    hdr = read_videomeg_vid(filename);
    checkUniqueLabels = false;
    
  case 'video'
    hdr = read_video(filename);
    checkUniqueLabels = false;
    
  otherwise
    if exist(headerformat, 'file')
      % attempt to run "headerformat" as a function, this allows the user to specify an external reading function
      % this is also used for bids_tsv, biopac_acq, motion_c3d, opensignals_txt, qualisys_tsv, sccn_xdf, and possibly others
      hdr = feval(headerformat, filename);
    elseif strcmp(fallback, 'biosig') && ft_hastoolbox('BIOSIG', 1)
      try
        % there is no guarantee that biosig can read it
        hdr = read_biosig_header(filename);
      catch
        ft_error('unsupported header format "%s"', headerformat);
      end
    else
      ft_error('unsupported header format "%s"', headerformat);
    end
    
end % switch headerformat


% Sometimes, the not all labels are correctly filled in by low-level reading functions. See for example bug #1572.
% First, make sure that there are enough (potentially empty) labels:
if numel(hdr.label) < hdr.nChans
  ft_warning('low-level reading function did not supply enough channel labels');
  hdr.label{hdr.nChans} = [];
end
% Now, replace all empty labels with new name:
if any(cellfun(@isempty, hdr.label))
  ft_warning('channel labels should not be empty, creating unique labels');
  hdr.label = fixlabels(hdr.label);
end

if checkUniqueLabels
  if length(hdr.label)~=length(unique(hdr.label))
    % all channels must have unique names
    ft_warning('all channels must have unique labels, creating unique labels');
    megflag = ft_chantype(hdr, 'meg');
    eegflag = ft_chantype(hdr, 'eeg');
    for i=1:hdr.nChans
      sel = find(strcmp(hdr.label{i}, hdr.label));
      if length(sel)>1
        % renaming the first instance is particularly disruptive when the channels are
        % part of standard MEG or EEG channel set, so that should be avoided
        if any(megflag(sel))
          sel = setdiff(sel, sel(find(megflag(sel), 1)));
        elseif any(eegflag(sel))
          sel = setdiff(sel, sel(find(eegflag(sel), 1)));
        else
          sel = sel(2:end);
        end
        for j=1:length(sel)
          % add a number to the original channel name
          hdr.label{sel(j)} = sprintf('%s-%d', hdr.label{sel(j)}, j);
        end
      end
    end
  end
end

% as of November 2011, the header is supposed to include the channel type (see FT_CHANTYPE,
% e.g. meggrad, megref, eeg) and the units of each channel (see FT_CHANUNIT, e.g. uV, fT)

if ~isfield(hdr, 'chantype') && checkUniqueLabels
  % use a helper function which has some built in intelligence
  hdr.chantype = ft_chantype(hdr);
end % for

if ~isfield(hdr, 'chanunit') && checkUniqueLabels
  % use a helper function which has some built in intelligence
  hdr.chanunit = ft_chanunit(hdr);
end % for

% ensure that the output grad is according to the latest definition
if isfield(hdr, 'grad')
  hdr.grad = ft_datatype_sens(hdr.grad);
end

% ensure that the output elec is according to the latest definition
if isfield(hdr, 'elec')
  hdr.elec = ft_datatype_sens(hdr.elec);
end

% ensure that the output opto is according to the latest definition
if isfield(hdr, 'opto')
  try
    hdr.opto = ft_datatype_sens(hdr.opto);
  catch
    % the NIRS optode structure is incomplete when reading/converting it from Homer files
    ft_warning('optode structure is not compliant with FT_DATATYPE_SENS');
  end
end

if (strcmp(readbids, 'yes') || strcmp(readbids, 'ifmakessense')) && isbids
  % the BIDS sidecar files overrule the information that is present in the file header itself
  try
    if exist('data_json', 'var')
      hdr.Fs = data_json.SamplingFrequency;
    end
    if exist('channels_tsv', 'var')
      assert(length(channels_tsv.name)  == hdr.nChans, 'number of channels is not consistent with the BIDS channels.tsv');
      assert(length(channels_tsv.type)  == hdr.nChans, 'number of channels is not consistent with the BIDS channels.tsv');
      assert(length(channels_tsv.units) == hdr.nChans, 'number of channels is not consistent with the BIDS channels.tsv');
      hdr.label     = channels_tsv.name;
      hdr.chantype  = channels_tsv.type;
      hdr.chanunit  = channels_tsv.units;
    end
    if exist('electrodes_tsv', 'var')
      hdr.elec         = [];
      hdr.elec.label   = electrodes_tsv.name;
      hdr.elec.elecpos = [electrodes_tsv.x electrodes_tsv.y electrodes_tsv.z];
    end
    if exist('optodes_tsv', 'var')
      hdr.opto         = [];
      hdr.opto.label   = optodes_tsv.name;
      hdr.opto.optopos = [optodes_tsv.x optodes_tsv.y optodes_tsv.z];
    end
  catch ME
    if strcmp(readbids, 'yes')
      ft_error(ME.message);
    else
      ft_warning(ME.message);
    end
  end % catch errors
end % if readbids and isbids

% ensure that these are column arrays and that they do not have empty entries
hdr.label = fixlabels(hdr.label);
if isfield(hdr, 'chantype'), hdr.chantype = fixchantype(hdr.chantype); end
if isfield(hdr, 'chanunit'), hdr.chanunit = fixchanunit(hdr.chanunit); end

% ensure that these are double precision and not integers, otherwise
% subsequent computations that depend on these might be messed up
hdr.Fs          = double(hdr.Fs);
hdr.nSamples    = double(hdr.nSamples);
hdr.nSamplesPre = double(hdr.nSamplesPre);
hdr.nTrials     = double(hdr.nTrials);
hdr.nChans      = double(hdr.nChans);

if inflated
  % compressed file has been unzipped on the fly, clean up
  if strcmp(headerformat, 'brainvision_vhdr')
    % don't delete the header file yet, ft_read_data might still need it
    % the files will be cleaned up by ft_read_data
  else
    delete(filename);
  end
end

if cache && exist(headerfile, 'file')
  % put the header in the cache
  cacheheader = hdr;
  % update the header details (including time stamp, size and name)
  cacheheader.details = dir(headerfile);
  % fprintf('added header to cache\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the file size in bytes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [siz] = filesize(filename)
l = dir(filename);
if l.isdir
  ft_error('"%s" is not a file', filename);
end
siz = l.bytes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the file size in bytes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr] = recursive_read_header(filename)
[p, f, x] = fileparts(filename);
ls = dir(filename);
ls = ls(~strcmp({ls.name}, '.'));  % exclude this directory
ls = ls(~strcmp({ls.name}, '..')); % exclude parent directory
for i=1:length(ls)
  % make sure that the directory listing includes the complete path
  ls(i).name = fullfile(filename, ls(i).name);
end
lst = {ls.name};
hdr = cell(size(lst));
sel = zeros(size(lst));
for i=1:length(lst)
  % read the header of each individual file
  try
    thishdr = ft_read_header(lst{i});
    if isstruct(thishdr)
      thishdr.filename = lst{i};
    end
  catch
    thishdr = [];
    ft_warning(lasterr);
    fprintf('while reading %s\n\n', lst{i});
  end
  if ~isempty(thishdr)
    hdr{i} = thishdr;
    sel(i) = true;
  else
    sel(i) = false;
  end
end
sel = logical(sel(:));
hdr = hdr(sel);
tmp = {};
for i=1:length(hdr)
  if isstruct(hdr{i})
    tmp = cat(1, tmp, hdr(i));
  elseif iscell(hdr{i})
    tmp = cat(1, tmp, hdr{i}{:});
  end
end
hdr = tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to fill in empty labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = fixlabels(labels)
for i = find(cellfun(@isempty, {labels{:}}))
  labels{i} = sprintf('%d', i);
end
labels = labels(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to fill in empty chantype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = fixchantype(labels)
sel = cellfun(@isempty, labels);
labels(sel) = {'unknown'};
labels = labels(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to fill in empty chanunit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = fixchanunit(labels)
sel = cellfun(@isempty, labels);
labels(sel) = {'unknown'};
labels = labels(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this is shared with DATA2BIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tsv = read_tsv(filename)
tsv = readtable(filename, 'Delimiter', 'tab', 'FileType', 'text', 'TreatAsEmpty', 'n/a', 'ReadVariableNames', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this is shared with DATA2BIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function json = read_json(filename)
ft_hastoolbox('jsonlab', 1);
json = loadjson(filename);
json = ft_struct2char(json); % convert strings into char-arrays
