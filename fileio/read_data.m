function [dat] = read_data(filename, varargin);

% READ_DATA reads electrophysiological data from a variety of EEG,
% MEG and LFP files and represents it in a common data-independent
% format. The supported formats are listed in the accompanying
% READ_HEADER function.
%
% Use as
%   dat = read_data(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'header'         header structure, see READ_HEADER
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
% See also READ_HEADER, READ_EVENT, WRITE_DATA, WRITE_EVENT

% Copyright (C) 2003-2007, Robert Oostenveld, F.C. Donders Centre
%
% $Log: read_data.m,v $
% Revision 1.89  2009/10/16 07:31:18  roboos
% renamed chieti into itab for consistency with other formats
%
% Revision 1.88  2009/10/13 10:12:51  roboos
% added support for chieti_raw
%
% Revision 1.87  2009/10/08 11:13:40  roevdmei
% added support for nmc_archive_k
%
% Revision 1.86  2009/10/07 12:45:32  roevdmei
% minor spelling correction
%
% Revision 1.85  2009/04/29 10:52:17  jansch
% incorporated handling of unsegmented egi simple binaries
%
% Revision 1.84  2009/03/23 12:09:07  vlalit
% Minor changes to make the ctf_old option fully functional.
%
% Revision 1.83  2009/03/13 07:12:09  roboos
% also use read_brainvision_eeg low-level function for *.seg files, obsoleted the read_brainvision_seg function
%
% Revision 1.82  2009/03/02 10:44:38  roboos
% switched default for fif files to use the MNE reading routines in case of neuromag_fif
% the user can make his own choise by specifying the format as neuromag_mne (for the MNE routines) or neuromag_mex (for the meg-pd mex files)
%
% Revision 1.81  2009/02/25 12:59:30  marvger
% removed calls to filetype and replaced with strcmp(dataformat,x) for
% efficiency
%
% Revision 1.80  2009/02/13 08:02:21  roboos
% ensure that the requested sample and trial numbers are integers
%
% Revision 1.79  2009/02/12 11:47:23  vlalit
% Added support for neuro prax (eldith) EEG format based on functions from the manufacturer
%  used with permission from the company's representative Mr. Klaus Schellhorn.
%
% Revision 1.78  2009/02/09 14:21:00  roboos
% added inport of micromed_trc data
%
% Revision 1.77  2009/02/09 13:35:16  roboos
% implemented efficient caching for bci2000,
% it should be initiated in read_header, subsequently read_data and read_event will reuse the details from the header
%
% Revision 1.76  2009/02/09 12:42:27  vlalit
% Added a check for discontinuous boundary in the case of neuromag_mne evoked data.
%
% Revision 1.75  2009/02/06 10:12:20  roboos
% incorporated the latest suggestions of Laurence Hunt for neuromag_mne
%
% Revision 1.74  2009/02/04 09:09:59  roboos
% fixed filename to headerfile/datafile cvonversion in case of ctf_old
%
% Revision 1.73  2009/01/23 16:18:36  roboos
% cleaned up neuromag_mne
%
% Revision 1.72  2009/01/23 10:32:55  vlalit
% New reader for Neuromag fif format using the MNE toolbox (http://www.nmr.mgh.harvard.edu/martinos/userInfo/data/sofMNE.php)  implemented by Laurence Hunt.
%
% Revision 1.71  2009/01/20 21:50:20  roboos
% use mxDeserialize instead of eval(string) in mysql
%
% Revision 1.70  2009/01/19 15:05:47  roboos
% added skeleton support for reading fif files using mne functions
%
% Revision 1.69  2009/01/06 09:12:03  roboos
% use true/false instead of 1/0
%
% Revision 1.68  2008/12/16 21:25:57  roboos
% removed the backward compatibility handling of the read_fcdc_data input arguments, these are now done in read_fcdc_data (i.e. keep it clean)
%
% Revision 1.67  2008/12/15 13:10:38  roboos
% fixed bug in biosig fallback, header would be read instead of data
%
% Revision 1.66  2008/12/01 14:49:59  roboos
% ensure that input arguments are double precision and not integers, otherwise the subsequent computations will be messed up (learned this in Lyon)
%
% Revision 1.65  2008/11/13 21:50:11  roboos
% also read 4d data in case ChannelUnitsPerBit is missing from header (give warning and set calibration to 1)
%
% Revision 1.64  2008/11/13 21:20:58  roboos
% use dataformat=ns_cnt16/32 to specify 16/32 bit format for reading neuroscan cnt files
%
% Revision 1.63  2008/11/02 10:59:41  roboos
% some more changes for ctf_ds in case of empty path
%
% Revision 1.62  2008/11/02 10:42:25  roboos
% improved handling of empty path in case of ctf dataset
%
% Revision 1.61  2008/09/29 21:46:02  roboos
% Implemented data caching in a data-format independent manner, using fetch_data and a persistent variable.
% Not yet suitable for inclusion in fileio release, hence the default is not to use caching.
%
% Revision 1.60  2008/09/29 08:37:44  roboos
% fixed bug when reading short segments of CTF data (errors were given on screen, so the bug was apparent)
%
% Revision 1.59  2008/09/25 11:53:48  roboos
% more efficient handling for CTF in real continuous datasets (i.e. one long trial) and in case all samples are within the same trial
% detected and fixed a bug that would case faulure of the code when reading a data segment that extends over more than two trials
%
% Revision 1.58  2008/09/24 16:26:17  roboos
% swiched from old fcdc import routines for CTF to the p-files supplied by CTF
% these new reading routines support synthetic gradients
% the format 'ctf_new' is not supported any more, because that is now the default
%
% Revision 1.57  2008/09/24 07:01:49  roboos
% fixed begsample for neuralynx_ncs (thanks to Martin)
%
% Revision 1.56  2008/07/24 08:44:20  roboos
% added initial support for nimh_cortex, not yet complete
%
% Revision 1.55  2008/07/01 16:23:02  roboos
% added read_combined_data (new implementation)
%
% Revision 1.54  2008/07/01 12:59:42  roboos
% explicit blockread 1 for ns_cnt
%
% Revision 1.53  2008/06/26 15:50:56  roboos
% tread ctf_new just as ctf_ds w.r.t. mapping of the filenames
%
% Revision 1.52  2008/06/20 07:25:56  roboos
% added check for presence of BCI2000 load_bcidat mex file
%
% Revision 1.51  2008/06/18 08:24:33  roboos
% added support for BCI2000
%
% Revision 1.50  2008/06/10 10:22:12  jansch
% removed sparse from calibration for 4d data, this is not necessary. added
% possibility for filnames including '.' for 4d data
%
% Revision 1.49  2008/05/29 13:54:52  roboos
% also work when no path is specified
%
% Revision 1.48  2008/05/29 13:51:11  roboos
% use strcmp instead of strmatch, thanks to Marinka
%
% Revision 1.47  2008/05/29 07:30:42  roboos
% small change in renaming header/data and filename in case of ctf, this prevents a warning if the res4 or meg4 are not positioned in a xxx.ds directory (see email Jo)
%
% Revision 1.46  2008/05/27 16:12:26  vlalit
% Changed type name to ced_spike6mat
%
% Revision 1.45  2008/05/27 11:58:20  vlalit
% Added support of Matlab files exported from Spike 6
%
% Revision 1.44  2008/05/15 15:10:56  roboos
% added ctf_new implementation, using p-files, this supports synthetic gradients
% some changes to the filename handling, merged nihm2grad into ctf2grad
%
% Revision 1.43  2008/05/02 17:45:10  vlalit
% Some bug fixes after testing the SPM fileo code
%
% Revision 1.41  2008/04/28 19:27:11  roboos
% fixed dimord and selection of eeglab set data
%
% Revision 1.40  2008/04/21 11:50:52  roboos
% added support for egi_sbin, thanks to Joseph Dien
%
% Revision 1.39  2008/04/18 14:07:45  roboos
% added eeglab_set
%
% Revision 1.38  2008/04/11 07:23:15  roboos
% updated docu
%
% Revision 1.37  2008/04/10 09:34:51  roboos
% added fallback option for biosig, implemented biosig also for edf
%
% Revision 1.36  2008/04/09 16:50:02  roboos
% added fallback option to biosig (not default)
%
% Revision 1.35  2008/04/09 14:10:34  roboos
% added placeholder for biosig (not yet implemented)
%
% Revision 1.34  2008/04/09 10:09:45  roboos
% pass the hdr.orig to the low-level brainvision readers
%
% Revision 1.33  2008/02/19 10:08:13  roboos
% added support for fcdc_buffer
%
% Revision 1.32  2008/01/31 20:12:51  roboos
% On line 322 the cell element 4D was added to the existing case
% designed for the handling of 4D_bti data.  This is necessary to
% allow this identified filetype to cause execution of this case.
% [thanks to Gavin]
%
% Revision 1.31  2008/01/10 12:57:34  roboos
% give explicit errors with msgid FILEIO:Something
%
% Revision 1.30  2007/12/17 13:03:52  roboos
% Vladimir found and fixed some bugs pertaining to the nexstim_nxe format
%
% Revision 1.29  2007/12/17 08:24:28  roboos
% added support for nexstim_nxe, thanks to Vladimir
% the low-level code has not been tested by myself
%
% Revision 1.28  2007/12/12 16:50:15  roboos
% added support for neuralynx_bin
%
% Revision 1.27  2007/11/07 10:49:06  roboos
% cleaned up the reading and writing from/to mysql database, using db_xxx helper functions (see mysql directory)
%
% Revision 1.26  2007/11/05 17:01:42  roboos
% added implementation for fcdc_mysql, not yet finished
%
% Revision 1.25  2007/09/24 15:14:02  roboos
% fixed bug in calibration of 4D/bti data, which is different for float than for short format data files (thanks to Nathan)
%
% Revision 1.24  2007/09/13 09:55:42  roboos
% use read_biosemi_bdf instead of openbdf/readbdf
%
% Revision 1.23  2007/08/01 12:24:40  roboos
% updated comments
%
% Revision 1.22  2007/08/01 09:57:01  roboos
% moved all code related to ctf shared memory to seperate functions
%
% Revision 1.21  2007/07/27 12:24:52  roboos
% implemented support for ctf_shm
% removed a double check for the presence of the file
%
% Revision 1.20  2007/07/19 14:49:21  roboos
% switched the default reader for nex files from read_nex_data to read_plexon_nex, the old one is still supported if explicitely mentioned as data/headerformat
% added support for multiple fragments of continuous AD data in nex files, holes are filled with NaNs
%
% Revision 1.19  2007/07/04 13:20:51  roboos
% added support for egi_egis/egia, thanks to Joseph Dien
%
% Revision 1.18  2007/07/03 16:10:48  roboos
% added a hack for 32 bit neuroscan format (hdr.nsdf=16|32), this should actually be be done using autodetection
%
% Revision 1.17  2007/07/03 15:53:46  roboos
% switched from using Cristian Wienbruchs BTi toolbox to a new ascii header reading function (read_bti_m4d)
%
% Revision 1.16  2007/06/13 08:06:21  roboos
% updated help
%
% Revision 1.15  2007/04/16 16:06:50  roboos
% fixed bug in selection of trials from an epoched file, when samples are requested (thanks to Vladimir)
%
% Revision 1.14  2007/03/26 12:41:18  roboos
% deal with continuous channels in nex files  that have different starting samples/timestamps
% small modification in plexon_plx to account fo the change in the API of the underlying function
%
% Revision 1.13  2007/03/21 17:24:01  roboos
% added plexon_ds
%
% Revision 1.12  2007/03/19 17:03:14  roboos
% added chanindx for read_neuralynx_dma
% implemented an alternative reader for nex files (read_plexon_nex), the old reader is still the default
%
% Revision 1.11  2007/02/21 09:53:44  roboos
% changed the rearrangement of the data in case of neuralynx_ncs to reflect the changed representation (2D numeric array instead of 1D cell array)
%
% Revision 1.10  2007/01/09 09:37:48  roboos
% added neuralynx_nte, added spike channels for plexon_plx
%
% Revision 1.9  2007/01/04 17:12:36  roboos
% implemented plexon_plx, only continuous channels
%
% Revision 1.8  2007/01/04 12:21:36  roboos
% updated the neuralynx section to reflect the new reading functions and to use read_neuralynx_cds
%
% Revision 1.7  2006/12/04 10:38:12  roboos
% added support for ns_avg
% fixed long-outstanding problem for reading multiple trials from ns_eeg data
%
% Revision 1.6  2006/09/18 21:47:54  roboos
% implemented support for fcdc_matbin, i.e. a dataset consisting of a matlab file with header and events and a seperate binary datafile
% added smart default for reading all data when no selection given
%
% Revision 1.5  2006/09/18 14:22:54  roboos
% implemented support for 4D-BTi dataformat
%
% Revision 1.4  2006/08/28 10:12:22  roboos
% use seperate filetype_check_extension instead of (missing) check_extension subfunction (thanks to Thilo for reporting this bug)
%
% Revision 1.3  2006/06/19 10:31:47  roboos
% translate continuous option into checkboundary (complements), added documentation
%
% Revision 1.2  2006/06/19 08:14:17  roboos
% updated documentation
%
% Revision 1.1  2006/06/07 09:32:20  roboos
% new implementation based on the read_fcdc_xxx functions, now with
% variable (key-val) input arguments, changed the control structure
% in the rpobram (switch instead of ifs), allow the user to specify
% the file format, allow the user to specify either a sample selection
% or a block selection. The reading functionality should not have
% changed compared to the read_fcdc_xxx versions.
%

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
  dataformat = filetype(filename);
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
if ~isempty(endsample) && mod(endsample, 1)
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
  case 'fcdc_matbin'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.mat']);
    datafile   = fullfile(path, [file '.bin']);
  case 'nmc_archive_k'
    [path, file, ext] = fileparts(filename);
    headerfile = [path '/' file 'newparams.txt'];
    if isempty(headerformat)
        headerformat = 'nmc_archive_k';
    end
    if isempty(hdr)
        hdr = read_header(headerfile, 'headerformat', headerformat);
    end
    datafile = filename;
  otherwise
    % convert filename into filenames, assume that the header and data are the same
    datafile   = filename;
    headerfile = filename;
end

if ~strcmp(filename, datafile) && ~ismember(dataformat, {'ctf_ds', 'ctf_old'})
  filename   = datafile;                % this function will read the data
  dataformat = filetype(filename);      % update the filetype
end

% for backward compatibility, default is to check when it is not continous
if isempty(checkboundary)
  checkboundary = ~keyval('continuous', varargin);
end

% read the header if it is not provided
if isempty(hdr)
  hdr = read_header(filename, 'headerformat', headerformat);
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
  if cache && isempty(cachedata)
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
      dat = fetch_data(cachedata, 'begsample', begsample', 'endsample', endsample);
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
    elseif length(chanindx)==hdr.nChans
      % all channels have been selected
      % nothing to do
    else
      % select the desired channel(s)
      dat = dat(chanindx,:);
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
    hastoolbox('BCI2000', 1);
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
    orig = read_besa_avr(filename);
    dat  = orig.data(chanindx, begsample:endsample);

  case 'besa_swf'
    % BESA source waveform
    orig = read_besa_swf(filename);
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

  case {'biosig', 'edf'}
    % use the biosig toolbox if available
    hastoolbox('BIOSIG', 1);
    dat = read_biosig_data(filename, hdr, begsample, endsample, chanindx);

  case {'brainvision_eeg', 'brainvision_dat', 'brainvision_seg'}
    dat = read_brainvision_eeg(filename, hdr.orig, begsample, endsample);
    dat = dat(chanindx,:);	% select the desired channels

  case 'ced_son'
    % chek the availability of the required low-level toolbox
    hastoolbox('neuroshare', 1);
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
    % check the presence of the required low-level toolbox
    hastoolbox('lc-libs', 1);
    chansel = hdr.orig.chansel;  % these are the channels that are visible to fieldtrip
    % read the data using the dll
    dat = lcReadData(chansel, begsample, endsample, filename);
    % take the subset of channels that is selected by the user
    dat = dat(chanindx, :);
    
  case  'combined_ds'
    dat = read_combined_ds(filename, hdr, begsample, endsample, chanindx);

  case {'ctf_ds', 'ctf_meg4', 'ctf_res4'}
    % check that the required low-level toolbox is available
    hastoolbox('ctf', 1);
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
    hastoolbox('eegsf', 1);
    % read it using the CTF importer from the NIH and Daren Weber
    tmp = ctf_read_meg4(filename, hdr.orig, chanindx, 'all', begtrial:endtrial);
    dat = cat(3, tmp.data{:});
    % the data is shaped in a 3-D array
    dimord = 'samples_chans_trials';

  case 'ctf_shm'
    % contact Robert Oostenveld if you are interested in real-time acquisition on the CTF system
    % read the data from shared memory
    [dat, dimord] = read_shm_data(hdr, chanindx, begtrial, endtrial);

  case 'eeglab_set'
    dat = read_eeglabdata(filename, 'header', hdr, 'begtrial', begtrial, 'endtrial', endtrial, 'chanindx', chanindx);
    dimord = 'chans_samples_trials';

  case 'spmeeg_mat'
    dat = read_spmeeg_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);

  case 'ced_spike6mat'
    dat = read_spike6mat_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);

  case 'eep_avr'
    % check that the required low-level toolbos ix available
    hastoolbox('eeprobe', 1);
    dat = read_eep_avr(filename);
    dat = dat.data(chanindx,begsample:endsample);       % select the desired channels and samples

  case 'eep_cnt'
    % check that the required low-level toolbos ix available
    hastoolbox('eeprobe', 1);
    dat = read_eep_cnt(filename, begsample, endsample);
    dat = dat.data(chanindx,:);                         % select the desired channels

  case 'fcdc_buffer'
    % read from a networked buffer for realtime analysis
    
    [host, port] = filetype_check_uri(filename);
    dat = buffer('get_dat', [begsample-1 endsample-1], host, port);  % indices should be zero-offset
    dat = dat.buf(chanindx,:);                                        % select the desired channels

  case 'fcdc_matbin'
    % multiplexed data in a *.bin file, accompanied by a matlab file containing the header
    offset        = begsample-1;
    numsamples    = endsample-begsample+1;
    samplesize    = 8;
    sampletype    = 'double';
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
    elseif length(chanindx)==hdr.nChans
      % all channels have been selected
      % nothing to do
    else
      % select the desired channel(s)
      dat = dat(chanindx,:);
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
    
  case 'micromed_trc'
    dat = read_micromed_trc(filename, begsample, endsample);
    dat = dat(chanindx,:);
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

  case 'nimh_cortex'
    keyboard

  case 'ns_avg'
    % NeuroScan average data
    orig = read_ns_avg(filename);
    dat  = orig.data(chanindx, begsample:endsample);

  case {'ns_cnt' 'ns_cnt16', 'ns_cnt32'}
    % Neuroscan continuous data
    sample1    = begsample-1;
    ldnsamples = endsample-begsample+1; % number of samples to read
    ldchan     = 1:hdr.nChans;          % must be row vector
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
    % read_ns_cnt originates from the EEGLAB package (loadcnt.m) but is
    % an old version since the new version is not compatible any more
    % all data is read, and only the relevant data is kept.
    if strcmp(dataformat, 'ns_cnt')
      tmp = read_ns_cnt(filename, 'sample1', sample1, 'ldnsamples', ldnsamples, 'ldchan', ldchan, 'blockread', 1);
    elseif strcmp(dataformat, 'ns_cnt16')
      tmp = read_ns_cnt(filename, 'sample1', sample1, 'ldnsamples', ldnsamples, 'ldchan', ldchan, 'blockread', 1, 'format', 16);
    elseif strcmp(dataformat, 'ns_cnt32')
      tmp = read_ns_cnt(filename, 'sample1', sample1, 'ldnsamples', ldnsamples, 'ldchan', ldchan, 'blockread', 1, 'format', 32);
    end
    dat = tmp.dat(chanoi,:);

  case 'ns_eeg'
    % Neuroscan epoched file
    tmp       = read_ns_eeg(filename, begtrial:endtrial);
    siz       = [(endtrial-begtrial+1) hdr.nChans hdr.nSamples];
    dat       = reshape(tmp.data, siz); % ensure 3-D array
    dat       = dat(:,chanindx,:);      % select channels
    dimord    = 'trials_chans_samples'; % selection using begsample and endsample will be done later

  case {'neuromag_fif' 'neuromag_mne'}
    % check that the required low-level toolbox is available
    hastoolbox('mne', 1);
    if (hdr.orig.iscontinuous)
      dat = fiff_read_raw_segment(hdr.orig.raw,begsample+hdr.nSamplesPre-1,endsample+hdr.nSamplesPre-1,chanindx);
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
    hastoolbox('meg-pd', 1);
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
    hastoolbox('yokogawa', 1);
    dat = read_yokogawa_data(filename, hdr, begsample, endsample, chanindx);

  case 'nmc_archive_k'
    dat = read_nmc_archive_k_data(filename, hdr, begsample, endsample, chanindx);

    
  otherwise
    if strcmp(fallback, 'biosig') && hastoolbox('BIOSIG', 1)
      dat = read_biosig_data(filename, hdr, begsample, endsample, chanindx);
    else
      error('unsupported data format');
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
  dat = reshape(dat, [nchans nsamples ntrials]); % convert into a 3-D array

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
