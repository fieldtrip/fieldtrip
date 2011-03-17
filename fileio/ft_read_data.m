function [dat] = ft_read_data(filename, varargin)

% FT_READ_DATA reads electrophysiological data from a variety of EEG,
% MEG and LFP files and represents it in a common data-independent
% format. The supported formats are listed in the accompanying
% FT_READ_HEADER function.
%
% Use as
%   dat = ft_read_data(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'header'         header structure, see FT_READ_HEADER
%   'begsample'      first sample to read
%   'endsample'      last sample to read
%   'begtrial'       first trial to read, mutually exclusive with begsample+endsample
%   'endtrial'       last trial to read, mutually exclusive with begsample+endsample
%   'chanindx'       list with channel indices to read
%   'checkboundary'  boolean, whether to check for reading segments over a trial boundary
%   'cache'          boolean, whether to use caching for multiple reads
%   'dataformat'     string
%   'headerformat'   string
%   'fallback'       can be empty or 'biosig' (default = [])
%
% This function returns a 2-D matrix of size Nchans*Nsamples for
% continuous data when begevent and endevent are specified, or a 3-D
% matrix of size Nchans*Nsamples*Ntrials for epoched or trial-based
% data when begtrial and endtrial are specified.
%
% See also FT_READ_HEADER, FT_READ_EVENT, FT_WRITE_DATA, FT_WRITE_EVENT

% Copyright (C) 2003-2010 Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

persistent cachedata     % for caching
persistent db_blob       % for fcdc_mysql

if isempty(db_blob)
  db_blob = 0;
end

% get the optional input arguments
hdr           = keyval('header',        varargin);
begsample     = keyval('begsample',     varargin);
endsample     = keyval('endsample',     varargin);
begtrial      = keyval('begtrial',      varargin);
endtrial      = keyval('endtrial',      varargin);
chanindx      = keyval('chanindx',      varargin);
checkboundary = keyval('checkboundary', varargin);
dataformat    = keyval('dataformat',    varargin);
headerformat  = keyval('headerformat',  varargin);
fallback      = keyval('fallback',      varargin);
cache         = keyval('cache',         varargin); if isempty(cache), cache = 0; end

% determine the filetype
if isempty(dataformat)
  dataformat = ft_filetype(filename);
end

% test whether the file or directory exists
if ~exist(filename, 'file') && ~strcmp(dataformat, 'ctf_shm') && ~strcmp(dataformat, 'fcdc_mysql') && ...
    ~strcmp(dataformat, 'fcdc_buffer')
  error('FILEIO:InvalidFileName', 'file or directory ''%s'' does not exist', filename);
end

% ensure that these are double precision and not integers, otherwise the subsequent computations will be messed up
begsample = double(begsample);
endsample = double(endsample);
begtrial  = double(begtrial);
endtrial  = double(endtrial);

% ensure that the requested sample and trial numbers are integers
if ~isempty(begsample) && mod(begsample, 1)
  warning('rounding "begsample" to the nearest integer');
  begsample = round(begsample);
end
if ~isempty(endsample) && ~isinf(endsample) && mod(endsample, 1)
  warning('rounding "endsample" to the nearest integer');
  endsample = round(endsample);
end
if ~isempty(begtrial) && mod(begtrial, 1)
  warning('rounding "begtrial" to the nearest integer');
  begtrial = round(begtrial);
end
if ~isempty(endtrial) && mod(endtrial, 1)
  warning('rounding "endtrial" to the nearest integer');
  endtrial = round(endtrial);
end

switch dataformat
  case '4d_pdf'
    datafile   = filename;
    headerfile = [datafile '.m4d'];
    sensorfile = [datafile '.xyz'];
  case {'4d_m4d', '4d_xyz'}
    datafile   = filename(1:(end-4)); % remove the extension
    headerfile = [datafile '.m4d'];
    sensorfile = [datafile '.xyz'];
  case '4d'
    [path, file, ext] = fileparts(filename);
    datafile   = fullfile(path, [file,ext]);
    headerfile = fullfile(path, [file,ext]);
    configfile = fullfile(path, 'config');
  case {'ctf_ds', 'ctf_old'}
    % convert CTF filename into filenames
    [path, file, ext] = fileparts(filename);
    if any(strcmp(ext, {'.res4' '.meg4', '.1_meg4' '.2_meg4' '.3_meg4' '.4_meg4' '.5_meg4' '.6_meg4' '.7_meg4' '.8_meg4' '.9_meg4'}))
      filename = path;
      [path, file, ext] = fileparts(filename);
    end
    if isempty(path) && isempty(file)
      % this means that the dataset was specified as the present working directory, i.e. only with '.'
      filename = pwd;
      [path, file, ext] = fileparts(filename);
    end
    headerfile = fullfile(filename, [file '.res4']);
    datafile   = fullfile(filename, [file '.meg4']);
    if length(path)>3 && strcmp(path(end-2:end), '.ds')
      filename = path; % this is the *.ds directory
    end
  case 'ctf_meg4'
    [path, file, ext] = fileparts(filename);
    if isempty(path)
      path = pwd;
    end
    headerfile = fullfile(path, [file '.res4']);
    datafile   = fullfile(path, [file '.meg4']);
    if length(path)>3 && strcmp(path(end-2:end), '.ds')
      filename = path; % this is the *.ds directory
    end
  case 'ctf_res4'
    [path, file, ext] = fileparts(filename);
    if isempty(path)
      path = pwd;
    end
    headerfile = fullfile(path, [file '.res4']);
    datafile   = fullfile(path, [file '.meg4']);
    if length(path)>3 && strcmp(path(end-2:end), '.ds')
      filename = path; % this is the *.ds directory
    end
  case 'brainvision_vhdr'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    if exist(fullfile(path, [file '.eeg']))
      datafile   = fullfile(path, [file '.eeg']);
    elseif exist(fullfile(path, [file '.seg']))
      datafile   = fullfile(path, [file '.seg']);
    elseif exist(fullfile(path, [file '.dat']))
      datafile   = fullfile(path, [file '.dat']);
    end
  case 'brainvision_eeg'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    datafile   = fullfile(path, [file '.eeg']);
  case 'brainvision_seg'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    datafile   = fullfile(path, [file '.seg']);
  case 'brainvision_dat'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    datafile   = fullfile(path, [file '.dat']);
  case 'itab_raw'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.raw.mhd']);
    datafile   = fullfile(path, [file '.raw']);
  case 'fcdc_matbin'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.mat']);
    datafile   = fullfile(path, [file '.bin']);
  case 'fcdc_buffer_offline'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '/header']);
    datafile = fullfile(path, [file '/samples']);
  case {'tdt_tsq' 'tdt_tev'}
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.tsq']);
    datafile   = fullfile(path, [file '.tev']);
  case 'nmc_archive_k'
    [path, file, ext] = fileparts(filename);
    headerfile = [path '/' file 'newparams.txt'];
    if isempty(headerformat)
      headerformat = 'nmc_archive_k';
    end
    if isempty(hdr)
      hdr = ft_read_header(headerfile, 'headerformat', headerformat);
    end
    datafile = filename;
  otherwise
    % convert filename into filenames, assume that the header and data are the same
    datafile   = filename;
    headerfile = filename;
end

if ~strcmp(filename, datafile) && ~ismember(dataformat, {'ctf_ds', 'ctf_old', 'fcdc_buffer_offline'})
  filename   = datafile;                % this function will read the data
  dataformat = ft_filetype(filename);      % update the filetype
end

% for backward compatibility, default is to check when it is not continous
if isempty(checkboundary)
  checkboundary = ~keyval('continuous', varargin);
end

% read the header if it is not provided
if isempty(hdr)
  hdr = ft_read_header(filename, 'headerformat', headerformat);
end

% set the default channel selection, which is all channels
if isempty(chanindx)
  chanindx = 1:hdr.nChans;
end

% read untill the end of the file if the endsample is "inf"
if any(isinf(endsample)) && any(endsample>0)
  endsample = hdr.nSamples*hdr.nTrials;
end

% test whether the requested channels can be accomodated
if min(chanindx)<1 || max(chanindx)>hdr.nChans
  error('FILEIO:InvalidChanIndx', 'selected channels are not present in the data');
end

% test whether the requested data segment is not outside the file
if any(begsample<1)
  error('FILEIO:InvalidBegSample', 'cannot read data before the begin of the file');
elseif any(endsample>(hdr.nSamples*hdr.nTrials))
  error('FILEIO:InvalidEndSample', 'cannot read data after the end of the file');
end

requesttrials  = isempty(begsample) && isempty(endsample);
requestsamples = isempty(begtrial)  && isempty(endtrial);

if cache && requesttrials
  error('caching is not supported when reading trials')
end

if isempty(begsample) && isempty(endsample) && isempty(begtrial) && isempty(endtrial)
  % neither samples nor trials are specified, set the defaults to read the complete data trial-wise (also works for continuous)
  requestsamples = 0;
  requesttrials  = 1;
  begtrial       = 1;
  endtrial       = hdr.nTrials;
end

% set the default, which is to assume that it is should check for boundaries when samples are requested
if isempty(checkboundary) && requesttrials
  checkboundary = false;
elseif isempty(checkboundary) && requestsamples
  checkboundary = true;
end

if requesttrials && requestsamples
  error('you cannot select both trials and samples at the same time');
elseif requesttrials
  % this allows for support using a continuous reader
  if isinf(hdr.nSamples) && begtrial==1
    begsample = 1;                             % computing it here does not work (0*inf=nan)
  else
    begsample = (begtrial-1)*hdr.nSamples + 1; % compute it the normal way
  end
  endsample = (endtrial  )*hdr.nSamples;
elseif requestsamples
  % this allows for support using a trial-based reader
  begtrial = floor((begsample-1)/hdr.nSamples)+1;
  endtrial = floor((endsample-1)/hdr.nSamples)+1;
else
  error('you should either specify begin/end trial or begin/end sample');
end

% test whether the requested data segment does not extend over a discontinuous trial boundary
if checkboundary && hdr.nTrials>1
  if begtrial~=endtrial
    error('requested data segment extends over a discontinuous trial boundary');
  end
end

if strcmp(dataformat, 'bci2000_dat')
  % caching for BCI2000 is handled in the main section and in read_header
else
  % implement the caching in a data-format independent way
  if cache && (isempty(cachedata) || ~isequal(cachedata.label,hdr.label(chanindx)))
    % create a new FieldTrip raw data structure that will hold the data
    cachedata.label = hdr.label(chanindx);
    cachedata.fsample = hdr.Fs;
    cachedata.time    = {};
    cachedata.trial   = {};
    cachedata.cfg     = [];
    cachedata.cfg.trl = zeros(0,3);
  elseif cache && ~isempty(cachedata)
    % try to fetch the requested segment from the cache
    try
      dat = ft_fetch_data(cachedata, 'begsample', begsample', 'endsample', endsample);
      % fprintf('caching succeeded\n');
      return
    catch
      % fprintf('caching failed\n');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the data with the low-level reading function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch dataformat

  case {'4d' '4d_pdf', '4d_m4d', '4d_xyz'}
    [fid,message] = fopen(datafile,'rb','ieee-be');
    % determine the type and size of the samples
    sampletype = lower(hdr.orig.Format);
    switch sampletype
      case 'short'
        samplesize = 2;
      case 'long'
        samplesize = 4;
      case 'float'
        samplesize = 4;
      case 'double'
        samplesize = 8;
      otherwise
        error('unsupported data format');
    end
    % 4D/BTi MEG data is multiplexed, can be epoched/discontinuous
    offset     = (begsample-1)*samplesize*hdr.nChans;
    numsamples = endsample-begsample+1;
    gain       = hdr.orig.ChannelGain;
    if isfield(hdr.orig, 'ChannelUnitsPerBit')
      upb = hdr.orig.ChannelUnitsPerBit;
    else
      warning('cannot determine ChannelUnitsPerBit');
      upb = 1;
    end
    % jump to the desired data
    fseek(fid, offset, 'cof');
    % read the desired data
    if length(chanindx)==1
      % read only one channel
      fseek(fid, (chanindx-1)*samplesize, 'cof');                                  % seek to begin of channel
      dat = fread(fid, numsamples, ['1*' sampletype], (hdr.nChans-1)*samplesize)'; % read one channel, skip the rest
    else
      % read all channels
      dat = fread(fid, [hdr.nChans, numsamples], sampletype);
    end
    fclose(fid);
    if length(chanindx)==1
      % only one channel was selected, which is managed by the code above
      % nothing to do
    elseif ~isequal(chanindx(:)', 1:hdr.nChans)
      dat = dat(chanindx,:);  % select the desired channels
    else
      % all channels have been selected
      % nothing to do
    end
    % determine how to calibrate the data
    switch sampletype
      case {'short', 'long'}
        % include both the gain values and the integer-to-double conversion in the calibration
        calib = diag(gain(chanindx) .* upb(chanindx));
      case {'float', 'double'}
        % only include the gain values in the calibration
        calib = diag(gain(chanindx));
      otherwise
        error('unsupported data format');
    end
    % calibrate the data
    dat = calib*dat;

  case 'bci2000_dat'
    % this requires the load_bcidat mex file to be present on the path
    ft_hastoolbox('BCI2000', 1);
    % this is inefficient, since it reads the complete data
    if isfield(hdr.orig, 'signal') && isfield(hdr.orig, 'states')
      % assume that the complete data is stored in the header, this speeds up subsequent read operations
      signal        = hdr.orig.signal;
      states        = hdr.orig.states;
      parameters    = hdr.orig.parameters;
      total_samples = hdr.orig.total_samples;
    else
      [signal, states, parameters, total_samples] = load_bcidat(filename);
    end
    % apply the callibration from AD units to uV
    dat = double(signal(begsample:endsample,chanindx)');
    for i=chanindx(:)'
      dat(i,:) = dat(i,:).* parameters.SourceChGain.NumericValue(i) + parameters.SourceChOffset.NumericValue(i);
    end
    dimord = 'chans_samples';

  case 'besa_avr'
    % BESA average data
    orig = readBESAavr(filename);
    dat  = orig.Data(chanindx, begsample:endsample);

  case 'besa_swf'
    % BESA source waveform
    orig = readBESAswf(filename);
    dat  = orig.data(chanindx, begsample:endsample);

  case {'biosemi_bdf', 'bham_bdf'}
    % this uses a mex file for reading the 24 bit data
    dat = read_biosemi_bdf(filename, hdr, begsample, endsample, chanindx);

  case {'biosemi_old'}
    % this uses the openbdf and readbdf functions that I copied from the EEGLAB toolbox
    epochlength = hdr.orig.Head.SampleRate(1);
    % it has already been checked in read_header that all channels have the same sampling rate
    begepoch = floor((begsample-1)/epochlength) + 1;
    endepoch = floor((endsample-1)/epochlength) + 1;
    nepochs  = endepoch - begepoch + 1;
    orig = openbdf(filename);
    dat  = zeros(length(chanindx),nepochs*epochlength);
    for i=begepoch:endepoch
      % read and concatenate all required data epochs
      [orig, buf] = readbdf(orig, i, 0);
      if size(buf,2)~=hdr.nChans || size(buf,1)~=epochlength
        error('error reading selected data from bdf-file');
      else
        dat(:,((i-begepoch)*epochlength+1):((i-begepoch+1)*epochlength)) = buf(:,chanindx)';
      end
    end
    begsample = begsample - (begepoch-1)*epochlength;  % correct for the number of bytes that were skipped
    endsample = endsample - (begepoch-1)*epochlength;  % correct for the number of bytes that were skipped
    dat = dat(:, begsample:endsample);
    % close the file between seperate read operations
    fclose(orig.Head.FILE.FID);

  case {'biosig', 'gdf'}
    % use the biosig toolbox if available
    ft_hastoolbox('BIOSIG', 1);
    dat = read_biosig_data(filename, hdr, begsample, endsample, chanindx);

  case {'brainvision_eeg', 'brainvision_dat', 'brainvision_seg'}
    dat = read_brainvision_eeg(filename, hdr.orig, begsample, endsample);
    if ~isequal(chanindx(:)', 1:hdr.nChans)
      dat = dat(chanindx,:);  % select the desired channels
    end

  case 'ced_son'
    % chek the availability of the required low-level toolbox
    ft_hastoolbox('neuroshare', 1);
    % use the reading function supplied by Gijs van Elswijk
    %
    % CED ADC is done in sequence, thus the analog channels
    % do not share the same time axis. This is _ignored_ here.
    %
    % Set the READ_CED_SON parameter 'readtimestamps' to
    % 'yes' to get time axis for each data channel returned.
    % This time information can be used to do
    % a temporal realignment of the data.
    tmp = read_ced_son(filename,'readdata','yes',...
      'begsample',begsample,...
      'endsample',endsample,...
      'channels',chanindx);
    dat = cell2mat(tmp.data');

  case  'itab_raw'
    if ismember(hdr.orig.data_type, [0 1 2])
      % big endian
      fid = fopen(datafile, 'rb', 'ieee-be');
    elseif ismember(hdr.orig.data_type, [3 4 5])
      % little endian
      fid = fopen(datafile, 'rb', 'ieee-le');
    else
      error('unsuppported data_type in itab format');
    end

    % skip the ascii header
    fseek(fid, hdr.orig.start_data, 'bof');

    if ismember(hdr.orig.data_type, [0 3])
      % short
      fseek(fid, (begsample-1)*hdr.orig.nchan*2, 'cof');
      dat = fread(fid, [hdr.orig.nchan endsample-begsample+1], 'int16');
    elseif ismember(hdr.orig.data_type, [1 4])
      % long
      fseek(fid, (begsample-1)*hdr.orig.nchan*4, 'cof');
      dat = fread(fid, [hdr.orig.nchan endsample-begsample+1], 'int32');
    elseif ismember(hdr.orig.data_type, [2 5])
      % float
      fseek(fid, (begsample-1)*hdr.orig.nchan*4, 'cof');
      dat = fread(fid, [hdr.orig.nchan endsample-begsample+1], 'float');
    else
      error('unsuppported data_type in itab format');
    end
    % these are the channels that are visible to fieldtrip
    chansel = 1:hdr.orig.nchan;  
    tmp = [hdr.orig.ch(chansel).calib];
    tmp = tmp(:);
    tmp(tmp==0) = 1;
    dat = dat ./ tmp(:,ones(1,size(dat,2)));
    % select the subset of visible channels that the user requested
    if ~isequal(chanindx(:)', 1:hdr.nChans)
      dat = dat(chanindx,:);  % select the desired channels
    end

  case  'combined_ds'
    dat = read_combined_ds(filename, hdr, begsample, endsample, chanindx);

  case {'ctf_ds', 'ctf_meg4', 'ctf_res4'}
    % check that the required low-level toolbox is available
    ft_hastoolbox('ctf', 1);
    % this returns SQUIDs in T, EEGs in V, ADC's and DACS in V, HLC channels in m, clock channels in s.
    if begtrial==endtrial
      % specify selection as 3x1 vector
      trlbegsample = begsample - hdr.nSamples*(begtrial-1); % within the trial
      trlendsample = endsample - hdr.nSamples*(begtrial-1); % within the trial
      dat = getCTFdata(hdr.orig, [trlbegsample; trlendsample; begtrial], chanindx, 'T', 'double');
      dimord = 'samples_chans';
    else
      % specify selection as 1xN vector
      dat = getCTFdata(hdr.orig, [begtrial:endtrial], chanindx, 'T', 'double');
      dimord = 'samples_chans_trials';
    end

  case  {'ctf_old', 'read_ctf_meg4'}
    % read it using the open-source matlab code that originates from CTF and that was modified by the FCDC
    dat = read_ctf_meg4(datafile, hdr.orig, begsample, endsample, chanindx);

  case 'ctf_read_meg4'
    % check that the required low-level toolbox is available
    ft_hastoolbox('eegsf', 1);
    % read it using the CTF importer from the NIH and Darren Weber
    tmp = ctf_read_meg4(filename, hdr.orig, chanindx, 'all', begtrial:endtrial);
    dat = cat(3, tmp.data{:});
    % the data is shaped in a 3-D array
    dimord = 'samples_chans_trials';

  case 'ctf_shm'
    % contact Robert Oostenveld if you are interested in real-time acquisition on the CTF system
    % read the data from shared memory
    [dat, dimord] = read_shm_data(hdr, chanindx, begtrial, endtrial);

  case 'dataq_wdq'
    dat = read_wdq_data(filename, hdr.orig, begsample, endsample, chanindx);
    
  case 'eeglab_set'
    dat = read_eeglabdata(filename, 'header', hdr, 'begtrial', begtrial, 'endtrial', endtrial, 'chanindx', chanindx);
    dimord = 'chans_samples_trials';

  case 'spmeeg_mat'
    dat = read_spmeeg_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);

  case 'ced_spike6mat'
    dat = read_spike6mat_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);

  case {'edf'}
    % this reader is largely similar to the one for bdf
    % it uses a mex file for reading the 16 bit data
    dat = read_edf(filename, hdr, begsample, endsample, chanindx);

  case 'eep_avr'
    % check that the required low-level toolbos ix available
    ft_hastoolbox('eeprobe', 1);
    dat = read_eep_avr(filename);
    dat = dat.data(chanindx,begsample:endsample);       % select the desired channels and samples

  case 'eep_cnt'
    % check that the required low-level toolbos ix available
    ft_hastoolbox('eeprobe', 1);
    dat = read_eep_cnt(filename, begsample, endsample);
    dat = dat.data(chanindx,:);                         % select the desired channels

  case 'eyelink_asc'
    if isfield(hdr.orig, 'dat')
      % this is inefficient, since it keeps the complete data in memory
      % but it does speed up subsequent read operations without the user
      % having to care about it
      asc = hdr.orig;
    else
      asc = read_eyelink_asc(filename);
    end
    dat = asc.dat(chanindx,begsample:endsample);

  case 'fcdc_buffer'
    % read from a networked buffer for realtime analysis
    [host, port] = filetype_check_uri(filename);
    dat = buffer('get_dat', [begsample-1 endsample-1], host, port);  % indices should be zero-offset
    dat = dat.buf(chanindx,:);                                        % select the desired channels

  case 'fcdc_buffer_offline'
    % read from a offline FieldTrip buffer data files
    dat = read_buffer_offline_data(datafile, hdr, [begsample endsample]);
    if ~isequal(chanindx(:)', 1:hdr.nChans)
      dat = dat(chanindx,:);  % select the desired channels
    end    

  case 'fcdc_matbin'
    % multiplexed data in a *.bin file, accompanied by a matlab file containing the header
    offset        = begsample-1;
    numsamples    = endsample-begsample+1;
    if isfield(hdr, 'precision'),
      sampletype    = hdr.precision;
    else
      sampletype    = 'double'; %original format without precision info in hdr is always in double
    end
    if strcmp(sampletype, 'single')
      samplesize    = 4;
    elseif strcmp(sampletype, 'double')
      samplesize    = 8;
    end
    [fid,message] = fopen(datafile,'rb','ieee-le');
    % jump to the desired data
    fseek(fid, offset*samplesize*hdr.nChans, 'cof');
    % read the desired data
    if length(chanindx)==1
      % read only one channel
      fseek(fid, (chanindx-1)*samplesize, 'cof');                                  % seek to begin of channel
      dat = fread(fid, numsamples, ['1*' sampletype], (hdr.nChans-1)*samplesize)'; % read one channel, skip the rest
    else
      % read all channels
      dat = fread(fid, [hdr.nChans, numsamples], sampletype);
    end
    fclose(fid);
    if length(chanindx)==1
      % only one channel was selected, which is managed by the code above
      % nothing to do
    elseif ~isequal(chanindx(:)', 1:hdr.nChans)
      dat = dat(chanindx,:);  % select the desired channels
    else
      % all channels have been selected
      % nothing to do
    end

  case 'fcdc_mysql'
    % read from a MySQL server listening somewhere else on the network
    db_open(filename);
    if db_blob
      error('not implemented');
    else
      for i=begtrial:endtrial
        s = db_select('fieldtrip.data', {'nChans', 'nSamples', 'data'}, i);
        dum{i-begtrial+1} = mxDeserialize(s.data);
      end
      dat = zeros(length(chanindx), s.nSamples, endtrial-begtrial+1);
      for i=begtrial:endtrial
        dat(:,:,i-begtrial+1) = dum{i-begtrial+1}(chanindx,:);
      end
      dimord = 'chans_samples_trials';
    end

  case {'egi_egia', 'egi_egis'}
    dat = read_egis_data(filename, hdr, begtrial, endtrial, chanindx);
    dimord = 'chans_samples_trials';

  case {'egi_sbin'}
    if mod(hdr.orig.header_array(1),2)==0,
      %unsegmented data contains only 1 trial, don't read the whole file
      dat = read_sbin_data(filename, hdr, begsample, endsample, chanindx);
      requestsamples = 0;
    else
      %segmented data
      dat = read_sbin_data(filename, hdr, begtrial, endtrial, chanindx);
    end
    dimord = 'chans_samples_trials';

  case {'egi_mff_bin'}
    % this is a file contained within a MFF package, which represents the complete dataset
    % better is to read the MFF package as a complete dataset instead of a single file
    blockhdr = hdr.orig;

    % the number of samples per block can be different
    % assume that all channels have the same sampling frequency and number of samples per block
    nsamples = zeros(size(blockhdr));
    for i=1:length(blockhdr)
      nsamples(i) = blockhdr(i).nsamples(1);
    end
    
    cumsamples = cumsum(nsamples);
    begblock = find(begsample<=cumsamples, 1, 'first');
    endblock = find(endsample<=cumsamples, 1, 'first');
    dat = read_mff_bin(filename, begblock, endblock);
    % select channels and concatenate in a matrix
    dat = cell2mat(dat(chanindx,:));

    % select the desired samples from the concatenated blocks
    if begblock==1
      prevsamples = 0;
    else
      prevsamples = cumsamples(begblock-1);
    end
    begsel = begsample-prevsamples;
    endsel = endsample-prevsamples;
    dat = dat(:,begsel:endsel);

  case 'micromed_trc'
    dat = read_micromed_trc(filename, begsample, endsample);
    if ~isequal(chanindx(:)', 1:hdr.nChans)
      dat = dat(chanindx,:);  % select the desired channels
    end
    dimord = 'chans_samples';

  case {'mpi_ds', 'mpi_dap'}
    [hdr, dat] = read_mpi_ds(filename);
    dat = dat(chanindx, begsample:endsample); % select the desired channels and samples

  case 'neuralynx_dma'
    dat = read_neuralynx_dma(filename, begsample, endsample, chanindx);

  case 'neuralynx_sdma'
    dat = read_neuralynx_sdma(filename, begsample, endsample, chanindx);

  case 'neuralynx_ncs'
    NRecords  = hdr.nSamples/512;
    begrecord = ceil(begsample/512);
    endrecord = ceil(endsample/512);
    % read the records that contain the desired samples
    ncs = read_neuralynx_ncs(filename, begrecord, endrecord);
    % cut out the desired samples
    begsample = begsample - (begrecord-1)*512;
    endsample = endsample - (begrecord-1)*512;
    % this also reshape the data from 512 X records into a linear array
    dat = ncs.dat(begsample:endsample);

  case 'neuralynx_nse'
    % read all records
    nse = read_neuralynx_nse(filename);
    % convert timestamps to samples
    sample = round((nse.TimeStamp - hdr.FirstTimeStamp)./hdr.TimeStampPerSample + 1);
    % select the timestamps that are between begin and endsample
    sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
    dat = zeros(1,endsample-begsample+1);
    dat(sample) = 1;

  case 'neuralynx_nte'
    % read all records
    nte = read_neuralynx_nte(filename);
    % convert timestamps to samples
    sample = round((nte.TimeStamp - hdr.FirstTimeStamp)./hdr.TimeStampPerSample + 1);
    % select the timestamps that are between begin and endsample
    sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
    dat = zeros(1,endsample-begsample+1);
    dat(sample) = 1;

  case {'neuralynx_ttl', 'neuralynx_tsl', 'neuralynx_tsh'}
    % single channel files
    dat = read_neuralynx_ttl(filename, begsample, endsample);

  case 'neuralynx_bin'
    % single channel files
    dat = read_neuralynx_bin(filename, begsample, endsample);

  case 'neuralynx_ds'
    dat = read_neuralynx_ds(filename, hdr, begsample, endsample, chanindx);

  case 'neuralynx_cds'
    dat = read_neuralynx_cds(filename, hdr, begsample, endsample, chanindx);

  case 'nexstim_nxe'
    dat = read_nexstim_nxe(filename, begsample, endsample, chanindx);

  case 'ns_avg'
    % NeuroScan average data
    orig = read_ns_avg(filename);
    dat  = orig.data(chanindx, begsample:endsample);

  case {'ns_cnt' 'ns_cnt16', 'ns_cnt32'}
    % Neuroscan continuous data
    sample1    = begsample-1;
    ldnsamples = endsample-begsample+1; % number of samples to read
    chanoi     = chanindx(:)';          % channels of interest
    if sample1<0
      error('begin sample cannot be for the beginning of the file');
    end
    % the hdr.nsdf was the initial fieldtrip hack to get 32 bit support, now it is realized using a extended dataformat string
    if     isfield(hdr, 'nsdf') && hdr.nsdf==16
      dataformat = 'ns_cnt16';
    elseif isfield(hdr, 'nsdf') && hdr.nsdf==32
      dataformat = 'ns_cnt32';
    end

    if strcmp(dataformat, 'ns_cnt')
      tmp = loadcnt(filename, 'sample1', sample1, 'ldnsamples', ldnsamples, 'blockread', 1);
    elseif strcmp(dataformat, 'ns_cnt16')
      tmp = loadcnt(filename, 'sample1', sample1, 'ldnsamples', ldnsamples, 'blockread', 1, 'dataformat', 'int16');
    elseif strcmp(dataformat, 'ns_cnt32')
      tmp = loadcnt(filename, 'sample1', sample1, 'ldnsamples', ldnsamples, 'blockread', 1, 'dataformat', 'int32');
    end
    dat = tmp.data(chanoi,:);

  case 'ns_eeg'
    % Neuroscan epoched file
    tmp       = read_ns_eeg(filename, begtrial:endtrial);
    siz       = [(endtrial-begtrial+1) hdr.nChans hdr.nSamples];
    dat       = reshape(tmp.data, siz); % ensure 3-D array
    dat       = dat(:,chanindx,:);      % select channels
    dimord    = 'trials_chans_samples'; % selection using begsample and endsample will be done later

  case {'neuromag_fif' 'neuromag_mne'}
    % check that the required low-level toolbox is available
    ft_hastoolbox('mne', 1);
    if (hdr.orig.iscontinuous)
      dat = fiff_read_raw_segment(hdr.orig.raw,begsample+hdr.orig.raw.first_samp-1,endsample+hdr.orig.raw.first_samp-1,chanindx);
      dimord = 'chans_samples';
    elseif (hdr.orig.isaverage)
      dat = cat(2, hdr.orig.evoked.epochs);            % concatenate all epochs, this works both when they are of constant or variable length
      if checkboundary
        trialnumber = [];
        for i = 1:numel(hdr.orig.evoked)
          trialnumber = [trialnumber i*ones(size(hdr.orig.evoked(i).times))];
        end
        if trialnumber(begsample) ~= trialnumber(endsample)
          error('requested data segment extends over a discontinuous trial boundary');
        end
      end
      dat = dat(chanindx, begsample:endsample);        % select the desired channels and samples
      dimord = 'chans_samples';
    elseif (hdr.orig.isepoched)
      error('Support for epoched *.fif data is not yet implemented.')
    end

  case 'neuromag_mex'
    % check that the required low-level toolbox is available
    ft_hastoolbox('meg-pd', 1);
    begtime = (begsample-1)/hdr.Fs;
    begepoch = floor((begsample-1)/hdr.nSamples) + 1;
    endepoch = floor((endsample-1)/hdr.nSamples) + 1;
    rawdata('any',filename);
    rawdata('goto', begtime);
    dat = [];
    for i=begepoch:endepoch
      [buf, status] = rawdata('next');
      if ~strcmp(status, 'ok')
        error('error reading selected data from fif-file');
      else
        dat(:,((i-begepoch)*hdr.nSamples+1):((i-begepoch+1)*hdr.nSamples)) = buf(chanindx,:);
      end
    end
    rawdata('close');
    begsample = begsample - (begepoch-1)*hdr.nSamples;  % correct for the number of bytes that were skipped
    endsample = endsample - (begepoch-1)*hdr.nSamples;  % correct for the number of bytes that were skipped
    dat = dat(:, begsample:endsample);

  case 'neuroprax_eeg'
    tmp = np_readdata(filename, hdr.orig, begsample - 1, endsample - begsample + 1, 'samples');
    dat = tmp.data';

  case 'plexon_ds'
    dat = read_plexon_ds(filename, hdr, begsample, endsample, chanindx);

  case 'plexon_ddt'
    dat = read_plexon_ddt(filename, begsample, endsample);
    dat = dat.data(chanindx,:);

  case {'read_nex_data'} % this is an alternative reader for nex files
    dat = read_nex_data(filename, hdr, begsample, endsample, chanindx);

  case {'read_plexon_nex' 'plexon_nex'} % this is the default reader for nex files
    dat = zeros(length(chanindx), endsample-begsample+1);
    for i=1:length(chanindx)
      if hdr.orig.VarHeader(chanindx(i)).Type==5
        % this is a continuous channel
        if hdr.orig.VarHeader(chanindx(i)).Count==1
          [nex, chanhdr] = read_plexon_nex(filename, 'header', hdr.orig, 'channel', chanindx(i), 'tsonly', 1);
          % the AD channel contains a single fragment
          % determine the sample offset into this fragment
          offset     = round(double(nex.ts-hdr.FirstTimeStamp)./hdr.TimeStampPerSample);
          chanbegsmp = begsample - offset;
          chanendsmp = endsample - offset;
          if chanbegsmp<1
            % the first sample of this channel is later than the beginning of the dataset
            % and we are trying to read the beginning of the dataset
            [nex, chanhdr] = read_plexon_nex(filename, 'header', hdr.orig, 'channel', chanindx(i), 'tsonly', 0, 'begsample', 1, 'endsample', chanendsmp);
            % padd the beginning of this channel with NaNs
            nex.dat = [nan(1,offset) nex.dat];
          else
            [nex, chanhdr] = read_plexon_nex(filename, 'header', hdr.orig, 'channel', chanindx(i), 'tsonly', 0, 'begsample', chanbegsmp, 'endsample', chanendsmp);
          end
          % copy the desired samples into the output matrix
          dat(i,:) = nex.dat;
        else
          % the AD channel contains multiple fragments
          [nex, chanhdr] = read_plexon_nex(filename, 'header', hdr.orig, 'channel', chanindx(i), 'tsonly', 0);
          % reconstruct the full AD timecourse with NaNs at all missing samples
          offset     = round(double(nex.ts-hdr.FirstTimeStamp)./hdr.TimeStampPerSample); % of each fragment, in AD samples
          nsample    = diff([nex.indx length(nex.dat)]);                                 % of each fragment, in AD samples
          % allocate memory to hold the complete continuous record
          cnt = nan(1, offset(end)+nsample(end));
          for j=1:length(offset)
            cntbegsmp  = offset(j)   + 1;
            cntendsmp  = offset(j)   + nsample(j);
            fragbegsmp = nex.indx(j) + 1;
            fragendsmp = nex.indx(j) + nsample(j);
            cnt(cntbegsmp:cntendsmp) = nex.dat(fragbegsmp:fragendsmp);
          end
          % copy the desired samples into the output matrix
          dat(i,:) = cnt(begsample:endsample);
        end
      elseif any(hdr.orig.VarHeader(chanindx(i)).Type==[0 1 3])
        % it is a neuron(0), event(1) or waveform(3) channel and therefore it has timestamps
        [nex, chanhdr] = read_plexon_nex(filename, 'header', hdr.orig, 'channel', chanindx(i), 'tsonly', 1);
        % convert the timestamps to samples
        sample = round(double(nex.ts - hdr.FirstTimeStamp)./hdr.TimeStampPerSample) + 1;
        % select only timestamps that are between begin and endsample
        sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
        for j=sample(:)'
          dat(i,j) = dat(i,j) + 1;
        end
      end
    end
    if any(isnan(dat(:)))
      warning('data has been padded with NaNs');
    end

  case 'plexon_plx'
    % determine the continuous channels
    contlabel = {hdr.orig.SlowChannelHeader.Name};
    for i=1:length(contlabel)
      contlabel{i} = deblank(contlabel{i});
    end
    [contindx, contsel]  = match_str(contlabel, hdr.label(chanindx));

    % determine the channels with spike waveforms
    spikelabel = {hdr.orig.ChannelHeader.Name};
    for i=1:length(spikelabel)
      spikelabel{i} = deblank(spikelabel{i});
    end
    [spikeindx, spikesel] = match_str(spikelabel, hdr.label(chanindx));

    if (length(contindx)+length(spikeindx))<length(chanindx)
      error('not all selected channels could be located in the data');
    end

    % allocate memory to hold all data
    dat = zeros(length(chanindx), endsample-begsample+1);

    % this is inefficient, since it reads all samples from each continuous channel
    % FIXME different continuous channels may start at a different timestamp
    for i=1:length(contsel)
      cont = read_plexon_plx(filename, 'header', hdr.orig, 'SlowChannelIndex', contindx(i), 'feedback', 1);
      dat(contsel(i),:) = cont(begsample:endsample);
    end

    % the timestamps of the spikes are in the header and do not have to be read
    for i=1:length(spikesel)
      % determine the timstamps of this channel
      sel = ([hdr.orig.DataBlockHeader.Type]==1 & [hdr.orig.DataBlockHeader.Channel]==hdr.orig.ChannelHeader(spikeindx(i)).Channel);
      tsl = [hdr.orig.DataBlockHeader(sel).TimeStamp];
      tsh = [hdr.orig.DataBlockHeader(sel).UpperByteOf5ByteTimestamp];
      ts  = timestamp_plexon(tsl, tsh); % use helper function, this returns an uint64 array
      % convert timestamps to samples
      sample = round(double(ts - hdr.FirstTimeStamp)./hdr.TimeStampPerSample + 1);
      % select only timestamps that are between begin and endsample
      sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
      for j=sample(:)'
        dat(spikesel(i),j) = dat(spikesel(i),j) + 1;
      end
    end

  case {'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw'}
    % check that the required low-level toolbox is available
    ft_hastoolbox('yokogawa', 1);
    dat = read_yokogawa_data(filename, hdr, begsample, endsample, chanindx);

  case 'nmc_archive_k'
    dat = read_nmc_archive_k_data(filename, hdr, begsample, endsample, chanindx);

  case 'neuroshare' % NOTE: still under development
    % check that the required neuroshare toolbox is available
    ft_hastoolbox('neuroshare', 1);

    tmp = read_neuroshare(filename, 'readanalog', 'yes', 'chanindx', chanindx, 'begsample', begsample, 'endsample', endsample);
    dat = tmp.analog.data';

  otherwise
    if strcmp(fallback, 'biosig') && ft_hastoolbox('BIOSIG', 1)
      dat = read_biosig_data(filename, hdr, begsample, endsample, chanindx);
    else
      error('unsupported data format (%s)', dataformat);
    end

end

if ~exist('dimord', 'var')
  dimord = 'chans_samples';  % almost all low-level readers return the data as 2D array
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reshape the 2-D or 3-D matrix to a common order of the dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch dimord
  case {'chans_samples', 'chans_samples_trials'}
    % nothing to do
  case 'samples_chans'
    dat = permute(dat, [2 1]);
    dimord = 'chans_samples';
  case 'samples_chans_trials'
    dat = permute(dat, [2 1 3]);
    dimord = 'chans_samples_trials';
  case 'trials_samples_chans'
    dat = permute(dat, [3 2 1]);
    dimord = 'chans_samples_trials';
  case 'trials_chans_samples'
    dat = permute(dat, [2 3 1]);
    dimord = 'chans_samples_trials';
  otherwise
    error('unexpected dimord');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between 3-D trial based and 2-D continuous output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if requesttrials  && strcmp(dimord, 'chans_samples')
  % reformat the continuous representation into trials
  nchans   = size(dat,1);
  nsamples = hdr.nSamples;
  ntrials  = size(dat,2)/hdr.nSamples;
  % in case of nchans=1 and ntrials=1, the reshaping into a 3D matrix results in the following warning
  %   Warning: ND-array output is being reshaped to a sparse 2D matrix.
  %          This behavior will change in a future release of MATLAB.
  if ntrials>1
    dat = reshape(dat, [nchans nsamples ntrials]); % convert into a 3-D array
  end

elseif requestsamples && strcmp(dimord, 'chans_samples_trials')
  % reformat the trials into a continuous representation
  nchans   = size(dat,1);
  nsamples = size(dat,2);
  ntrials  = size(dat,3);
  dat = reshape(dat, [nchans nsamples*ntrials]); % convert into a 2-D array
  % determine the selection w.r.t. the data as it is on disk
  begselection = (begtrial-1)*hdr.nSamples + 1;
  endselection = (endtrial  )*hdr.nSamples;
  % determine the selection w.r.t. the data that has been read in
  begselection2 = begsample - begselection + 1;
  endselection2 = endsample - begselection + 1;
  dat = dat(:,begselection2:endselection2);
end

if strcmp(dataformat, 'bci2000_dat')
  % caching for BCI2000 is handled in the main section and in read_header
else
  % implement caching in a data independent way
  if cache && requestsamples
    % add the new segment to the cache
    % FIMXE the cache size should be limited
    cachedata.cfg.trl(end+1,:) = [begsample endsample 0];
    cachedata.trial{end+1} = dat;
    cachedata.time{end+1} = (1:size(dat,2))/cachedata.fsample;
  end
end
