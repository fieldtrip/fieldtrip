function [dat] = ft_read_data(filename, varargin)

% FT_READ_DATA reads electrophysiological data from a variety of EEG, MEG and LFP
% files and represents it in a common data-independent format. The supported formats
% are listed in the accompanying FT_READ_HEADER function.
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
%   'chanunit'       cell-array with strings, the desired unit of each channel
%   'checkboundary'  boolean, whether to check for reading segments over a trial boundary
%   'checkmaxfilter' boolean, whether to check that maxfilter has been correctly applied (default = true)
%   'cache'          boolean, whether to use caching for multiple reads
%   'dataformat'     string
%   'headerformat'   string
%   'fallback'       can be empty or 'biosig' (default = [])
%   'blocking'       wait for the selected number of events (default = 'no')
%   'timeout'        amount of time in seconds to wait when blocking (default = 5)
%
% This function returns a 2-D matrix of size Nchans*Nsamples for continuous
% data when begevent and endevent are specified, or a 3-D matrix of size
% Nchans*Nsamples*Ntrials for epoched or trial-based data when begtrial
% and endtrial are specified.
%
% The list of supported file formats can be found in FT_READ_HEADER.
%
% See also FT_READ_HEADER, FT_READ_EVENT, FT_WRITE_DATA, FT_WRITE_EVENT

% Copyright (C) 2003-2016 Robert Oostenveld
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

persistent cachedata     % for caching
persistent db_blob       % for fcdc_mysql

if isempty(db_blob)
  db_blob = 0;
end

if iscell(filename)
  ft_warning(sprintf('concatenating data from %d files', numel(filename)));
  % this only works if the data is indexed by means of samples, not trials
  assert(isempty(ft_getopt(varargin, 'begtrial')));
  assert(isempty(ft_getopt(varargin, 'endtrial')));
  % use recursion to read data from multiple files
  
  hdr = ft_getopt(varargin, 'header');
  if isempty(hdr) || ~isfield(hdr, 'orig') || ~iscell(hdr.orig)
    for i=1:numel(filename)
      % read the individual file headers
      hdr{i}  = ft_read_header(filename{i}, varargin{:});
    end
  else
    % use the individual file headers that were read previously
    hdr = hdr.orig;
  end
  nsmp = nan(size(filename));
  for i=1:numel(filename)
    nsmp(i) = hdr{i}.nSamples*hdr{i}.nTrials;
  end
  offset = [0 cumsum(nsmp(1:end-1))];
  
  dat       = cell(size(filename));
  begsample = ft_getopt(varargin, 'begsample', 1);
  endsample = ft_getopt(varargin, 'endsample', sum(nsmp));
  
  for i=1:numel(filename)
    thisbegsample = begsample - offset(i);
    thisendsample = endsample - offset(i);
    if thisbegsample<=nsmp(i) && thisendsample>=1
      varargin = ft_setopt(varargin, 'header', hdr{i});
      varargin = ft_setopt(varargin, 'begsample', max(thisbegsample,1));
      varargin = ft_setopt(varargin, 'endsample', min(thisendsample,nsmp(i)));
      dat{i} = ft_read_data(filename{i}, varargin{:});
    else
      dat{i} = [];
    end
  end
  
  % return the concatenated data
  dat = cat(2, dat{:});
  return
end

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

% get the optional input arguments
hdr             = ft_getopt(varargin, 'header');
begsample       = ft_getopt(varargin, 'begsample');
endsample       = ft_getopt(varargin, 'endsample');
begtrial        = ft_getopt(varargin, 'begtrial');
endtrial        = ft_getopt(varargin, 'endtrial');
chanindx        = ft_getopt(varargin, 'chanindx');
checkboundary   = ft_getopt(varargin, 'checkboundary');
checkmaxfilter  = ft_getopt(varargin, 'checkmaxfilter', 'yes'); % this is only passed as varargin to FT_READ_HEADER
headerformat    = ft_getopt(varargin, 'headerformat');
fallback        = ft_getopt(varargin, 'fallback');
cache           = ft_getopt(varargin, 'cache', false);
dataformat      = ft_getopt(varargin, 'dataformat');
chanunit        = ft_getopt(varargin, 'chanunit');
timestamp       = ft_getopt(varargin, 'timestamp');

% this allows blocking reads to avoid having to poll many times for online processing
blocking         = ft_getopt(varargin, 'blocking', false);  % true or false
timeout          = ft_getopt(varargin, 'timeout', 5);       % seconds

% convert from 'yes'/'no' into boolean
blocking = istrue(blocking);

if isempty(dataformat)
  % only do the autodetection if the format was not specified
  dataformat = ft_filetype(filename);
end

if iscell(dataformat)
  % this happens for datasets specified as cell-array for concatenation
  dataformat = dataformat{1};
end

% test whether the file or directory exists
if ~any(strcmp(dataformat, {'fcdc_buffer', 'ctf_shm', 'fcdc_mysql'})) && ~exist(filename, 'file')
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

if strcmp(dataformat, 'compressed')
  % the file is compressed, unzip on the fly
  inflated   = true;
  filename   = inflate_file(filename);
  dataformat = ft_filetype(filename);
else
  inflated   = false;
end

% ensure that the headerfile and datafile are defined, which are sometimes different than the name of the dataset
[filename, headerfile, datafile] = dataset2files(filename, dataformat);

if ~strcmp(filename, datafile) && ~any(strcmp(dataformat, {'ctf_ds', 'ctf_old', 'fcdc_buffer_offline'}))
  filename   = datafile;                % this function will read the data
  dataformat = ft_filetype(filename);   % update the filetype
end

% for backward compatibility, default is to check when it is not continous
if isempty(checkboundary)
  checkboundary = ~ft_getopt(varargin, 'continuous');
end

% read the header if it is not provided
if isempty(hdr)
  hdr = ft_read_header(filename, 'headerformat', headerformat, 'chanindx', chanindx, 'checkmaxfilter', checkmaxfilter);
  if isempty(chanindx)
    chanindx = 1:hdr.nChans;
  end
else
  % set the default channel selection, which is all channels
  if isempty(chanindx)
    chanindx = 1:hdr.nChans;
  end
  % test whether the requested channels can be accomodated
  if min(chanindx)<1 || max(chanindx)>hdr.nChans
    error('FILEIO:InvalidChanIndx', 'selected channels are not present in the data');
  end
end

% read until the end of the file if the endsample is "inf"
if any(isinf(endsample)) && any(endsample>0)
  endsample = hdr.nSamples*hdr.nTrials;
end

% test whether the requested data segment is not outside the file
if any(begsample<1)
  error('FILEIO:InvalidBegSample', 'cannot read data before the begin of the file');
elseif any(endsample>(hdr.nSamples*hdr.nTrials)) && ~blocking
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

if any(strcmp(dataformat, {'bci2000_dat', 'eyelink_asc', 'gtec_mat', 'mega_neurone'}))
  % caching for these formats is handled in the main section and in ft_read_header
else
  % implement the caching in a data-format independent way
  if cache && (isempty(cachedata) || ~isequal(cachedata.label,hdr.label(chanindx)))
    % create a new FieldTrip raw data structure that will hold the data
    cachedata.label = hdr.label(chanindx);
    cachedata.fsample = hdr.Fs;
    cachedata.time    = {};
    cachedata.trial   = {};
    cachedata.cfg     = [];
    cachedata.sampleinfo = zeros(0,2);
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
% please maintain this list in alphabetical order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch dataformat
  
  case {'4d' '4d_pdf', '4d_m4d', '4d_xyz'}
    [fid, message] = fopen(datafile,'rb','ieee-be');
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
    dat = double(full(sparse(calib)*dat));
    
  case 'AnyWave'
    dat = read_ah5_data(filename, hdr, begsample, endsample, chanindx);
    
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
    
  case 'besa_besa'
    dat = read_besa_besa(filename, hdr, begsample, endsample, chanindx);
    
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
    
  case 'biosemi_old'
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
    % close the file between separate read operations
    fclose(orig.Head.FILE.FID);
    
  case {'biosig'}
    % this requires the biosig toolbox
    ft_hastoolbox('BIOSIG', 1);
    dat = read_biosig_data(filename, hdr, begsample, endsample, chanindx);
    
    
  case 'blackrock_nsx'
    % use the NPMK toolbox for the file reading
    ft_hastoolbox('NPMK', 1);
    
    % ensure that the filename contains a full path specification,
    % otherwise the low-level function fails
    [p,f,e] = fileparts(filename);
    if ~isempty(p)
      % this is OK
    elseif isempty(p)
      filename = which(filename);
    end
    orig = openNSx(filename, 'duration', [begsample endsample]);
    keyboard
    
  case {'brainvision_eeg', 'brainvision_dat', 'brainvision_seg'}
    dat = read_brainvision_eeg(filename, hdr.orig, begsample, endsample, chanindx);
    
  case 'bucn_nirs'
    dat = read_bucn_nirsdata(filename, hdr, begsample, endsample, chanindx);
    
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
    
  case 'combined_ds'
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
    % read it using the open-source MATLAB code that originates from CTF and that was modified by the FCDC
    dat = read_ctf_meg4(datafile, hdr.orig, begsample, endsample, chanindx);
    
  case 'ctf_read_meg4'
    % check that the required low-level toolbox is available
    ft_hastoolbox('eegsf', 1);
    % read it using the CTF importer from the NIH and Darren Weber
    tmp = ctf_read_meg4(fileparts(datafile), hdr.orig, chanindx, 'all', begtrial:endtrial);
    dat = cat(3, tmp.data{:});
    % the data is shaped in a 3-D array
    dimord = 'samples_chans_trials';
    
  case 'ctf_shm'
    % contact Robert Oostenveld if you are interested in real-time acquisition on the CTF system
    % read the data from shared memory
    [dat, dimord] = read_shm_data(hdr, chanindx, begtrial, endtrial);
    
  case 'ced_spike6mat'
    dat = read_spike6mat_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);
    
  case {'deymed_ini' 'deymed_dat'}
    % the data is stored in a binary *.dat file
    if isempty(hdr)
      hdr.orig = [];
    end
    dat = read_deymed_dat(datafile, hdr.orig, begsample, endsample);
    dat = dat(chanindx, :);
    
  case 'dataq_wdq'
    dat = read_wdq_data(filename, hdr.orig, begsample, endsample, chanindx);
    
  case 'eeglab_set'
    dat = read_eeglabdata(filename, 'header', hdr, 'begtrial', begtrial, 'endtrial', endtrial, 'chanindx', chanindx);
    dimord = 'chans_samples_trials';
    
  case 'eeglab_erp'
    dat = read_erplabdata(filename, 'header', hdr, 'begtrial', begtrial, 'endtrial', endtrial, 'chanindx', chanindx);
    dimord = 'chans_samples_trials';
    
  case 'emotiv_mat'
    % This is a MATLAB *.mat file that is created using the Emotiv MATLAB
    % example code. It contains a 25xNsamples matrix and some other stuff.
    dat = hdr.orig.data_eeg';
    dat = dat(chanindx, begsample:endsample);
    
  case {'egi_egia', 'egi_egis'}
    dat = read_egis_data(filename, hdr, begtrial, endtrial, chanindx);
    dimord = 'chans_samples_trials';
    
  case 'egi_sbin'
    if (bitand(hdr.orig.header_array(1),1) == 0) && ~((hdr.orig.header_array(14)==0) && (hdr.orig.header_array(15) > 1)),
      %unsegmented data contains only 1 trial, don't read the whole file
      dat = read_sbin_data(filename, hdr, begsample, endsample, chanindx);
      requestsamples = 0;
    else
      %segmented data
      dat = read_sbin_data(filename, hdr, begtrial, endtrial, chanindx);
    end
    dimord = 'chans_samples_trials';
    
  case {'egi_mff_v1' 'egi_mff'} % this is currently the default
    
    % The following represents the code that was written by Ingrid, Robert
    % and Giovanni to get started with the EGI mff dataset format. It might
    % not support all details of the file formats.
    % An alternative implementation has been provided by EGI, this is
    % released as fieldtrip/external/egi_mff and referred further down in
    % this function as 'egi_mff_v2'.
    
    % check if requested data contains multiple epochs and not segmented. If so, give error
    if isfield(hdr.orig.xml,'epochs') && length(hdr.orig.xml.epochs) > 1
      if hdr.nTrials ==1
        data_in_epoch = zeros(1,length(hdr.orig.xml.epochs));
        for iEpoch = 1:length(hdr.orig.xml.epochs)
          begsamp_epoch = hdr.orig.epochdef(iEpoch,1);
          endsamp_epoch = hdr.orig.epochdef(iEpoch,2);
          data_in_epoch(iEpoch) = length(intersect(begsamp_epoch:endsamp_epoch,begsample:endsample));
        end
        if sum(data_in_epoch>1) > 1
          warning('The requested segment from %i to %i is spread out over multiple epochs with possibly discontinuous boundaries', begsample, endsample);
        end
      end
    end
    
    % read in data in different signals
    binfiles = dir(fullfile(filename, 'signal*.bin'));
    if isempty(binfiles)
      error('FieldTrip:read_mff_header:nobin', ['could not find any signal.bin in ' filename_mff ])
    end
    % determine which channels are in which signal
    for iSig = 1:length(hdr.orig.signal)
      if iSig == 1
        chan2sig_ind(1:hdr.orig.signal(iSig).blockhdr(1).nsignals(1)) = iSig;
      else
        chan2sig_ind(end+1:end+1+hdr.orig.signal(iSig).blockhdr(1).nsignals(1)) = iSig;
      end
    end
    for iSig = 1:length(hdr.orig.signal)
      % adjust chanindx to match with current signal
      [dum1, dum2, chanind_sig] = intersect(chanindx, find(chan2sig_ind==iSig));
      if isempty(chanind_sig)
        % no channels requested from current signal
      else
        blockhdr = hdr.orig.signal(iSig).blockhdr;
        signalname = binfiles(iSig).name;
        fullsignalname = fullfile(filename, signalname);
        
        % the number of samples per block can be different
        % assume that all channels have the same sampling frequency and number of samples per block
        nsamples = zeros(size(blockhdr));
        for i=1:length(blockhdr)
          nsamples(i) = blockhdr(i).nsamples(1);
        end
        
        cumsamples = cumsum(nsamples);
        begblock = find(begsample<=cumsamples, 1, 'first');
        endblock = find(endsample<=cumsamples, 1, 'first');
        datsig = read_mff_bin(fullsignalname, begblock, endblock, chanind_sig);
        
        % concatenate in a matrix
        if exist('dat', 'var')
          dat{length(dat)+1} = cell2mat(datsig(:,:));
        else
          dat{1} = cell2mat(datsig(:,:));
        end
        % select the desired samples from the concatenated blocks
        if begblock==1
          prevsamples = 0;
        else
          prevsamples = cumsamples(begblock-1);
        end
        begsel = begsample-prevsamples;
        endsel = endsample-prevsamples;
        dat{end} = dat{end}(:,begsel:endsel);
      end
    end
    % concat signals
    dat = cat(1,dat{:});
    
    if hdr.nTrials > 1
      dat2=zeros(hdr.nChans,hdr.nSamples,hdr.nTrials);
      for i=1:hdr.nTrials
        dat2(:,:,i)=dat(:,hdr.orig.epochdef(i,1):hdr.orig.epochdef(i,2));
      end;
      dat=dat2;
    end
    
  case 'egi_mff_v2'
    % ensure that the EGI_MFF toolbox is on the path
    ft_hastoolbox('egi_mff', 1);
    % ensure that the JVM is running and the jar file is on the path
    %%%%%%%%%%%%%%%%%%%%%%
    %workaround for MATLAB bug resulting in global variables being cleared
    globalTemp=cell(0);
    globalList=whos('global');
    varList=whos;
    for i=1:length(globalList)
      eval(['global ' globalList(i).name ';']);
      eval(['globalTemp{end+1}=' globalList(i).name ';']);
    end;
    %%%%%%%%%%%%%%%%%%%%%%
    
    mff_setup;
    
    %%%%%%%%%%%%%%%%%%%%%%
    %workaround for MATLAB bug resulting in global variables being cleared
    varNames={varList.name};
    for i=1:length(globalList)
      eval([globalList(i).name '=globalTemp{i};']);
      if ~any(strcmp(globalList(i).name,varNames)) %was global variable originally out of scope?
        eval(['clear ' globalList(i).name ';']); %clears link to global variable without affecting it
      end;
    end;
    clear globalTemp globalList varNames varList;
    %%%%%%%%%%%%%%%%%%%%%%
    
    if isunix && filename(1)~=filesep
      % add the full path to the dataset directory
      filename = fullfile(pwd, filename);
    elseif ispc && filename(2)~=':'
      % add the full path, including drive letter
      filename = fullfile(pwd, filename);
    end
    % pass the header along to speed it up, it will be read on the fly in case it is empty
    dat = read_mff_data(filename, 'sample', begsample, endsample, chanindx, hdr);
    
  case 'edf'
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
    
    if blocking
      nsamples  = endsample; % indices should be zero-offset
      nevents   = 0;         % disable waiting for events
      available = buffer_wait_dat([nsamples nevents timeout], host, port);
      if available.nsamples<nsamples
        error('buffer timed out while waiting for %d samples', nsamples);
      end
    end
    
    dat = buffer('get_dat', [begsample-1 endsample-1], host, port);  % indices should be zero-offset
    dat = dat.buf(chanindx,:);                                       % select the desired channels
    
  case 'fcdc_buffer_offline'
    % read from a offline FieldTrip buffer data files
    dat = read_buffer_offline_data(datafile, hdr, [begsample endsample]);
    if ~isequal(chanindx(:)', 1:hdr.nChans)
      dat = dat(chanindx,:);  % select the desired channels
    end
    
  case 'fcdc_matbin'
    % multiplexed data in a *.bin file, accompanied by a MATLAB file containing the header
    offset        = begsample-1;
    numsamples    = endsample-begsample+1;
    if isfield(hdr, 'precision'),
      sampletype  = hdr.precision;
    else
      sampletype  = 'double'; %original format without precision info in hdr is always in double
    end
    if strcmp(sampletype, 'single')
      samplesize  = 4;
    elseif strcmp(sampletype, 'double')
      samplesize  = 8;
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
    % check that the required low-level toolbox is available
    ft_hastoolbox('mysql', 1);
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
      fileset = {filename};
      count = 0;
      while exist(sprintf('%s_%d%s', fullfile(p, f), count+1, x), 'file')
        fileset{end+1} = sprintf('%s_%d%s', fullfile(p, f), count+1, x);
        count = count+1;
      end
      
      % determine which parts have to be read from which file
      nSamples = [hdr.orig.nSamples] .* [hdr.orig.nTrials];
      fileBegSample = [0 cumsum(nSamples(1:end-1))]+1;
      fileEndSample = cumsum(nSamples);
      
      dat = cell(1,length(fileset));
      for i=1:length(fileset)
        if begsample<=fileEndSample(i) && endsample>=fileBegSample(i)
          % read a piece of data from this file
          thisBegSample = begsample - fileBegSample(i) + 1;
          thisEndSample = endsample - fileBegSample(i) + 1;
          thisBegSample = max(1,           thisBegSample);
          thisEndSample = min(nSamples(i), thisEndSample);
          dat{i} = read_biosig_data(fileset{i}, hdr.orig(i), thisBegSample, thisEndSample, chanindx);
        else
          dat{i} = zeros(length(chanindx),0);
        end
      end
      % concatenate the data from the different files
      dat = cat(2, dat{:});
      
    else
      % there is only a single file
      dat = read_biosig_data(filename, hdr, begsample, endsample, chanindx);
    end
    
  case 'gtec_mat'
    if isfield(hdr, 'orig')
      % these are remembered in the hdr.orig field for fast reading of subsequent segments
      log   = hdr.orig.log;
      names = hdr.orig.names;
    else
      % this is a simple MATLAB format, it contains a log and a names variable
      tmp = load(headerfile);
      log   = tmp.log;
      names = tmp.names;
    end
    dat = log(chanindx, begsample:endsample);
    dimord = 'chans_samples';
    
  case 'itab_raw'
    if any(hdr.orig.data_type==[0 1 2])
      % big endian
      fid = fopen(datafile, 'rb', 'ieee-be');
    elseif any(hdr.orig.data_type==[3 4 5])
      % little endian
      fid = fopen(datafile, 'rb', 'ieee-le');
    else
      error('unsuppported data_type in itab format');
    end
    
    % skip the ascii header
    fseek(fid, hdr.orig.start_data, 'bof');
    
    if any(hdr.orig.data_type==[0 3])
      % short
      fseek(fid, (begsample-1)*hdr.orig.nchan*2, 'cof');
      dat = fread(fid, [hdr.orig.nchan endsample-begsample+1], 'int16');
    elseif any(hdr.orig.data_type==[1 4])
      % long
      fseek(fid, (begsample-1)*hdr.orig.nchan*4, 'cof');
      dat = fread(fid, [hdr.orig.nchan endsample-begsample+1], 'int32');
    elseif any(hdr.orig.data_type==[2 5])
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
    
  case 'jaga16'
    fid = fopen(filename, 'r');
    fseek(fid, hdr.orig.offset + (begtrial-1)*hdr.orig.packetsize, 'bof');
    buf = fread(fid, (endtrial-begtrial+1)*hdr.orig.packetsize/2, 'uint16');
    fclose(fid);
    % the packet is 1396 bytes with timestamp or 1388 without
    packet = jaga16_packet(buf, hdr.orig.packetsize==1396);
    % Our amplifier was rated as +/- 5mV input signal range, and we use 16
    % bit ADC.  However when we actually measured the signal range in our
    % device the input range can go as high as +/- 6 mV.  In this case our
    % bit resolution is about 0.2uV/bit. (instead of 0.16uV/bit)
    calib  = 0.2;
    dat    = calib * packet.dat;
    dimord = 'chans_samples';
    
  case {'manscan_mb2', 'manscan_mbi'}
    [p, f, x] = fileparts(filename);
    filename  = fullfile(p, [f, '.mb2']);
    trlind = [];
    if isfield(hdr.orig, 'epochs') && ~isempty(hdr.orig.epochs)
      for i = 1:numel(hdr.orig.epochs)
        trlind = [trlind i*ones(1, diff(hdr.orig.epochs(i).samples) + 1)];
      end
      if checkboundary && (trlind(begsample)~=trlind(endsample))
        error('requested data segment extends over a discontinuous trial boundary');
      end
    else
      trlind = ones(1, hdr.nSamples);
    end
    
    iEpoch = unique(trlind(begsample:endsample));
    sfid = fopen(filename, 'r');
    dat  = zeros(hdr.nChans, endsample - begsample + 1);
    for i = 1:length(iEpoch)
      dat(:, trlind(begsample:endsample) == iEpoch(i)) =...
        in_fread_manscan(hdr.orig, sfid, iEpoch(i), ...
        [sum(trlind==iEpoch(i) & (1:length(trlind))<begsample) ...
        sum(trlind==iEpoch(i) & (1:length(trlind))<=endsample)-1]);
    end
    dat = dat(chanindx, :);
    
  case 'mega_neurone'
    % this is fast but memory inefficient, since the header contains all data and events
    if isfield(hdr.orig, 'data')
      NEURONE = hdr.orig;
    else
      % ensure that this external toolbox is on the path
      ft_hastoolbox('neurone', 1);
      if filename(end)~=filesep
        % it should end with a slash
        filename = [filename filesep];
      end
      NEURONE = readneurone(filename);
    end
    dat = NEURONE.data(chanindx, begsample:endsample);
    dimord = 'chans_samples';
    
  case 'micromed_trc'
    dat = read_micromed_trc(filename, begsample, endsample);
    if ~isequal(chanindx(:)', 1:hdr.nChans)
      dat = dat(chanindx,:);  % select the desired channels
    end
    dimord = 'chans_samples';
    
  case {'mpi_ds', 'mpi_dap'}
    [hdr, dat] = read_mpi_ds(filename);
    dat = dat(chanindx, begsample:endsample); % select the desired channels and samples
    
  case 'neuroscope_bin'
    switch hdr.orig.nBits
      case 16
        precision = 'int16';
      case 32
        precision = 'int32';
      otherwise
        error('unknown precision');
    end
    dat     = LoadBinary(filename, 'frequency', hdr.Fs, 'offset', begsample-1, 'nRecords', endsample-begsample, 'nChannels', hdr.orig.nChannels, 'channels', chanindx, 'precision', precision).';
    scaling = hdr.orig.voltageRange/hdr.orig.amplification/(2^hdr.orig.nBits); % scale to S.I. units, i.e. V
    dat     = scaling.*dat;
    
  case 'netmeg'
    % the data is in the same NetCDF file as the header and is cached in the header structure
    dat = hdr.orig.Var.waveforms;
    dat = dat(:,:,chanindx);
    % at the bottom of ft_read_data there is a general handling of the permuting and reshaping
    dimord = 'trials_samples_chans';
    
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
    if istrue(timestamp)
      ncs.dat = cast(ncs.dat, class(ncs.TimeStamp));
      d = ncs.TimeStamp(2:end)-ncs.TimeStamp(1:end-1);
      medianTimestampPerBlock  = median(double(d)); % to avoid influence of the gaps
      TimestampPerSample       = medianTimestampPerBlock/512; % divide by known block size
      cls = class(ncs.TimeStamp);
      % replace the data with the timestamp of each sample
      for i=1:512
        ncs.dat(i,:) = ncs.TimeStamp + cast((i-1)*TimestampPerSample,cls);
      end
    end
    % this selects samples and also reshape the data from 512*Nrecords into a linear array (row)
    dat = ncs.dat(begsample:endsample);
    dat = dat(:)';
    
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
    ft_hastoolbox('eeglab', 1);
    % Neuroscan continuous data
    sample1    = begsample-1;
    ldnsamples = endsample-begsample+1; % number of samples to read
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
      tmp = loadcnt(filename, 'sample1', sample1, 'ldnsamples', ldnsamples); % let loadcnt figure it out
    elseif strcmp(dataformat, 'ns_cnt16')
      tmp = loadcnt(filename, 'sample1', sample1, 'ldnsamples', ldnsamples, 'dataformat', 'int16');
    elseif strcmp(dataformat, 'ns_cnt32')
      tmp = loadcnt(filename, 'sample1', sample1, 'ldnsamples', ldnsamples, 'dataformat', 'int32');
    end
    dat = tmp.data(chanindx,:);
    
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
    elseif (hdr.orig.isepoched)
      data = permute(hdr.orig.epochs.data, [2 3 1]);  % Chan X Sample X Trials
      if requesttrials
        dat = data(chanindx, :, begtrial:endtrial);
      else
        dat = data(chanindx, begsample:endsample);  % reading over boundaries
      end
    elseif (hdr.orig.isaverage)
      assert(isfield(hdr.orig, 'evoked'), '%s does not contain evoked data', filename);
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
    
  case {'neurosim_ds' 'neurosim_signals'}
    [hdr, dat] = read_neurosim_signals(filename);
    if endsample>size(dat,2)
      warning('Simulation was not completed, reading in part of the data')
      endsample=size(dat,2);
    end
    dat = dat(chanindx,begsample:endsample);
    
  case 'neurosim_evolution'
    [hdr, dat] = read_neurosim_evolution(filename);
    if endsample>size(dat,2)
      warning('Simulation was not completed, reading in part of the data')
      endsample=size(dat,2);
    end
    dat = dat(chanindx,begsample:endsample);
    
  case 'neurosim_spikes'
    warning('Reading Neurosim spikes as continuous data, for better memory efficiency use spike structure provided by ft_read_spike instead.');
    spike = ft_read_spike(filename);
    cfg          = [];
    cfg.trialdef.triallength = inf;
    cfg.trialfun = 'ft_trialfun_general';
    cfg.trlunit='samples'; %ft_trialfun_general gives us samples, not timestamps
    
    cfg.datafile=filename;
    cfg.hdr = ft_read_header(cfg.datafile);
    warning('off','FieldTrip:ft_read_event:unsupported_event_format')
    cfg = ft_definetrial(cfg);
    warning('on','FieldTrip:ft_read_event:unsupported_event_format')
    spiketrl = ft_spike_maketrials(cfg,spike);
    
    dat=ft_checkdata(spiketrl,'datatype', 'raw', 'fsample', spiketrl.hdr.Fs);
    dat=dat.trial{1};
    
  case 'nmc_archive_k'
    dat = read_nmc_archive_k_data(filename, hdr, begsample, endsample, chanindx);
    
  case 'neuroshare' % NOTE: still under development
    % check that the required neuroshare toolbox is available
    ft_hastoolbox('neuroshare', 1);
    tmp = read_neuroshare(filename, 'readanalog', 'yes', 'chanindx', chanindx, 'begsample', begsample, 'endsample', endsample);
    dat = tmp.analog.data';
    
  case 'neuroprax_eeg'
    tmp = np_readdata(filename, hdr.orig, begsample - 1, endsample - begsample + 1, 'samples');
    dat = tmp.data(:,chanindx)';
    
  case 'oxy3'
    dat = read_artinis_oxy3(filename, hdr, begsample, endsample, chanindx);
    
  case 'plexon_ds'
    dat = read_plexon_ds(filename, hdr, begsample, endsample, chanindx);
    
  case 'plexon_ddt'
    dat = read_plexon_ddt(filename, begsample, endsample);
    dat = dat.data(chanindx,:);
    
  case {'plexon_nex' 'read_plexon_nex'} % this is the default reader for nex files
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
    
  case 'read_nex_data' % this is an alternative reader for nex files
    dat = read_nex_data(filename, hdr, begsample, endsample, chanindx);
    
  case 'riff_wave'
    dat = wavread(filename, [begsample endsample])';
    dat = dat(chanindx,:);
    
  case 'spmeeg_mat'
    dat = read_spmeeg_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);
    
  case 'tmsi_poly5'
    blocksize = hdr.orig.header.SamplePeriodsPerBlock;
    begtrial = floor((begsample-1)/blocksize) + 1;
    endtrial = floor((endsample-1)/blocksize) + 1;
    dat = read_tmsi_poly5(filename, hdr.orig, begtrial, endtrial);
    offset = (begtrial-1)*blocksize;
    % select the desired samples and channels
    dat = dat(chanindx, (begsample-offset):(endsample-offset));
    
  case 'videomeg_aud'
    dat = read_videomeg_aud(filename, hdr, begsample, endsample);
    dat = dat(chanindx,:);
    
  case 'videomeg_vid'
    dat = read_videomeg_vid(filename, hdr, begsample, endsample);
    dat = dat(chanindx,:);
    
  case {'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw'}
    % the data can be read with three toolboxes: Yokogawa MEG Reader, Maryland sqdread,
    % or Yokogawa MEG160 (old inofficial toolbox)
    % newest toolbox takes precedence over others.
    
    if ft_hastoolbox('yokogawa_meg_reader', 3); %stay silent if it cannot be added
      dat = read_yokogawa_data_new(filename, hdr, begsample, endsample, chanindx);
    elseif ft_hastoolbox('sqdproject', 2) % give warning if it cannot be added
      % channels are counted 0-based, samples are counted 1-based
      [dat, info] = sqdread(filename, 'channels', chanindx-1, 'samples', [begsample endsample]);
      dat = dat';
    else
      ft_hastoolbox('yokogawa', 1); % error if it cannot be added
      dat = read_yokogawa_data(filename, hdr, begsample, endsample, chanindx);
    end
    
  otherwise
    if strcmp(fallback, 'biosig') && ft_hastoolbox('BIOSIG', 1)
      dat = read_biosig_data(filename, hdr, begsample, endsample, chanindx);
    else
      error('unsupported data format (%s)', dataformat);
    end
end % switch dataformat

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
% convert the channel data to the desired units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(chanunit)
  if length(chanunit)~=length(chanindx)
    error('the number of channel units is inconsistent with the number of channels');
  end
  
  % determine the scaling factor for each channel
  scaling = cellfun(@ft_scalingfactor, hdr.chanunit(chanindx(:)), chanunit(:));
  
  switch dimord
    case 'chans_samples'
      for i=1:length(scaling)
        dat(i,:) = scaling(i) .* dat(i,:);
      end
    case'chans_samples_trials';
      for i=1:length(scaling)
        dat(i,:,:) = scaling(i) .* dat(i,:,:);
      end
    otherwise
      error('unexpected dimord');
  end % switch
end % if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between 3-D trial based and 2-D continuous output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if requesttrials  && strcmp(dimord, 'chans_samples')
  % reformat the continuous representation into trials
  nchans   = size(dat,1);
  nsamples = hdr.nSamples;
  ntrials  = size(dat,2)/hdr.nSamples;
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

if inflated
  % compressed file has been unzipped on the fly, clean up
  if strcmp(dataformat, 'brainvision_eeg')
    % delete the complete directory, including the header and marker file
    delete(fileparts(filename));
  else
    delete(filename);
  end
end

if strcmp(dataformat, 'bci2000_dat') || strcmp(dataformat, 'eyelink_asc') || strcmp(dataformat, 'gtec_mat')
  % caching for these formats is handled in the main section and in read_header
else
  % implement caching in a data independent way
  if cache && requestsamples
    % add the new segment to the cache
    % FIMXE the cache size should be limited
    cachedata.sampleinfo(end+1,:) = [begsample endsample];
    cachedata.trial{end+1} = dat;
    cachedata.time{end+1} = (1:size(dat,2))/cachedata.fsample;
  end
end
