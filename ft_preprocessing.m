function [data] = ft_preprocessing(cfg, data)

% FT_PREPROCESSING reads MEG and/or EEG data according to user-specified trials
% and applies several user-specified preprocessing steps to the signals.
%
% Use as
%   [data] = ft_preprocessing(cfg)
% or
%   [data] = ft_preprocessing(cfg, data)
%
% The first input argument "cfg" is the configuration structure, which contains all
% details for the dataset filename, trials and the preprocessing options.
%
% If you are calling FT_PREPROCESSING with only the configuration as first input
% argument and the data still has to be read from file, you should specify
%   cfg.dataset      = string with the filename
%   cfg.trl          = Nx3 matrix with the trial definition, see FT_DEFINETRIAL
%   cfg.padding      = length (in seconds) to which the trials are padded for filtering (default = 0)
%   cfg.padtype      = string, type of padding (default: 'data' padding or
%                      'mirror', depending on feasibility)
%   cfg.continuous   = 'yes' or 'no' whether the file contains continuous data
%                      (default is determined automatic)
%
% Instead of specifying the dataset in the configuration, you can also explicitly
% specify the name of the file containing the header information and the name of the
% file containing the data, using
%   cfg.datafile     = string with the filename
%   cfg.headerfile   = string with the filename
%
% If you are calling FT_PREPROCESSING with the second input argument "data", then
% that should contain data that was already read from file in a previous call to
% FT_PREPROCESSING. In that case only the configuration options below apply.
%
% The channels that will be read and/or preprocessed are specified with
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.chantype     = string or Nx1 cell-array with channel types to be read (only for NeuroOmega)
%
% The preprocessing options for the selected channels are specified with
%   cfg.lpfilter      = 'no' or 'yes'  lowpass filter (default = 'no')
%   cfg.hpfilter      = 'no' or 'yes'  highpass filter (default = 'no')
%   cfg.bpfilter      = 'no' or 'yes'  bandpass filter (default = 'no')
%   cfg.bsfilter      = 'no' or 'yes'  bandstop filter (default = 'no')
%   cfg.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform (default = 'no')
%   cfg.medianfilter  = 'no' or 'yes'  jump preserving median filter (default = 'no')
%   cfg.lpfreq        = lowpass  frequency in Hz
%   cfg.hpfreq        = highpass frequency in Hz
%   cfg.bpfreq        = bandpass frequency range, specified as [lowFreq highFreq] in Hz
%   cfg.bsfreq        = bandstop frequency range, specified as [low high] in Hz (or as Nx2 matrix for notch filter)
%   cfg.dftfreq       = line noise frequencies in Hz for DFT filter (default = [50 100 150])
%   cfg.lpfiltord     = lowpass  filter order (default set in low-level function)
%   cfg.hpfiltord     = highpass filter order (default set in low-level function)
%   cfg.bpfiltord     = bandpass filter order (default set in low-level function)
%   cfg.bsfiltord     = bandstop filter order (default set in low-level function)
%   cfg.lpfilttype    = digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')
%   cfg.hpfilttype    = digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')
%   cfg.bpfilttype    = digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')
%   cfg.bsfilttype    = digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')
%   cfg.lpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.hpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.bpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.bsfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.lpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.hpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.bpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.bsinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.lpfiltdf      = lowpass transition width (firws, overrides order, default set in low-level function)
%   cfg.hpfiltdf      = highpass transition width (firws, overrides order, default set in low-level function)
%   cfg.bpfiltdf      = bandpass transition width (firws, overrides order, default set in low-level function)
%   cfg.bsfiltdf      = bandstop transition width (firws, overrides order, default set in low-level function)
%   cfg.lpfiltwintype = lowpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.hpfiltwintype = highpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.bpfiltwintype = bandpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.bsfiltwintype = bandstop window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.lpfiltdev     = lowpass max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.hpfiltdev     = highpass max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.bpfiltdev     = bandpass max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.bsfiltdev     = bandstop max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.dftreplace    = 'zero' or 'neighbour', method used to reduce line noise, 'zero' implies DFT filter, 'neighbour' implies spectrum interpolation (default = 'zero')
%   cfg.dftbandwidth  = bandwidth of line noise frequencies, applies to spectrum interpolation, in Hz (default = [1 2 3])
%   cfg.dftneighbourwidth = bandwidth of frequencies neighbouring line noise frequencies, applies to spectrum interpolation, in Hz (default = [2 2 2])
%   cfg.plotfiltresp  = 'no' or 'yes', plot filter responses (firws, default = 'no')
%   cfg.usefftfilt    = 'no' or 'yes', use fftfilt instead of filter (firws, default = 'no')
%   cfg.medianfiltord = length of median filter (default = 9)
%   cfg.demean        = 'no' or 'yes', whether to apply baseline correction (default = 'no')
%   cfg.baselinewindow = [begin end] in seconds, the default is the complete trial (default = 'all')
%   cfg.detrend       = 'no' or 'yes', remove linear trend from the data (done per trial) (default = 'no')
%   cfg.polyremoval   = 'no' or 'yes', remove higher order trend from the data (done per trial) (default = 'no')
%   cfg.polyorder     = polynome order for poly trend removal (default = 2; note that all lower-order trends will also be removed when using cfg.polyremoval)
%   cfg.derivative    = 'no' or 'yes', computes the first order derivative of the data (default = 'no')
%   cfg.hilbert       = 'no', 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag' or 'angle' (default = 'no')
%   cfg.rectify       = 'no' or 'yes' (default = 'no')
%   cfg.precision     = 'single' or 'double' (default = 'double')
%   cfg.absdiff       = 'no' or 'yes', computes absolute derivative (i.e. first derivative then rectify)
%
% Preprocessing options that only apply to MEG data are
%   cfg.coordsys      = string, 'head' or 'dewar' (default = 'head')
%   cfg.coilaccuracy  = can be empty or a number (0, 1 or 2) to specify the accuracy (default = [])
%   cfg.coildeffile   = can be empty or a string to a custom coil_def.dat file (default = [])
%
% Preprocessing options that you should only use for EEG data are
%   cfg.reref         = 'no' or 'yes' (default = 'no')
%   cfg.refchannel    = cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
%   cfg.refmethod     = 'avg', 'median', 'rest', 'bipolar' or 'laplace' (default = 'avg')
%   cfg.groupchans    = 'yes' or 'no', should channels be rereferenced in separate groups for bipolar and laplace methods,
%                       this requires channnels to be named using an alphanumeric code, where letters represent the group
%                       and numbers represent the order of the channel whithin its group (default = 'no')
%   cfg.leadfield     = leadfield structure, this is required when cfg.refmethod='rest', see FT_PREPARE_LEADFIELD
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%   cfg.montage       = 'no' or a montage structure, see FT_APPLY_MONTAGE (default = 'no')
%
% Preprocessing options that you should only use when you are calling FT_PREPROCESSING with
% also the second input argument "data" are
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%
% Preprocessing options that you should only use when you are calling
% FT_PREPROCESSING with a single cfg input argument are
%   cfg.method        = 'trial' or 'channel', read data per trial or per channel (default = 'trial')
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_DEFINETRIAL, FT_REDEFINETRIAL, FT_APPENDDATA, FT_APPENDSPIKE

% Undocumented local options:
%   cfg.paddir     = direction of padding, 'left'/'right'/'both' (default = 'both')
%   cfg.artfctdef =
%   cfg.removemcg =
%   cfg.montage   = (in combination with meg-data in the input) applies montage to both data and grad-structure)
%
% You can use this function to read data from one format, filter it, and
% write it to disk in another format. The reading is done either as one
% long continuous segment or in multiple trials. This is achieved by
%   cfg.export.dataset    = string with the output file name
%   cfg.export.dataformat = string describing the output file format, see FT_WRITE_DATA

% Copyright (C) 2003-2024, Robert Oostenveld, SMI, FCDC
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% return immediately after distributed execution
if ~isempty(ft_getopt(cfg, 'distribute'))
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels', 'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamed',    {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',    {'blcwindow', 'baselinewindow'});
cfg = ft_checkconfig(cfg, 'renamed',    {'output', 'export'});

% set the defaults
cfg.method         = ft_getopt(cfg, 'method', 'trial');
cfg.channel        = ft_getopt(cfg, 'channel', 'all');
cfg.removemcg      = ft_getopt(cfg, 'removemcg', 'no');
cfg.removeeog      = ft_getopt(cfg, 'removeeog', 'no');
cfg.precision      = ft_getopt(cfg, 'precision', 'double');
cfg.padding        = ft_getopt(cfg, 'padding', 0);          % padding is only done when filtering
cfg.paddir         = ft_getopt(cfg, 'paddir', 'both');
cfg.montage        = ft_getopt(cfg, 'montage', 'no');
cfg.updatesens     = ft_getopt(cfg, 'updatesens', 'no');    % in case a montage or rereferencing is specified
cfg.dataformat     = ft_getopt(cfg, 'dataformat');          % is passed to low-level function, empty implies autodetection

% these options relate to the actual preprocessing, it is necessary to specify here because of padding
cfg.dftfilter      = ft_getopt(cfg, 'dftfilter', 'no');
cfg.lpfilter       = ft_getopt(cfg, 'lpfilter', 'no');
cfg.hpfilter       = ft_getopt(cfg, 'hpfilter', 'no');
cfg.bpfilter       = ft_getopt(cfg, 'bpfilter', 'no');
cfg.bsfilter       = ft_getopt(cfg, 'bsfilter', 'no');
cfg.medianfilter   = ft_getopt(cfg, 'medianfilter', 'no');
cfg.padtype        = ft_getopt(cfg, 'padtype', 'data');

% these options relate to the actual preprocessing, it is necessary to specify here because of channel selection
cfg.reref          = ft_getopt(cfg, 'reref', 'no');
cfg.refchannel     = ft_getopt(cfg, 'refchannel', {});
cfg.refmethod      = ft_getopt(cfg, 'refmethod', 'avg');
cfg.groupchans     = ft_getopt(cfg, 'groupchans', 'no');
cfg.implicitref    = ft_getopt(cfg, 'implicitref');

% construct the low-level options as key-value pairs, these are passed to FT_READ_HEADER and FT_READ_DATA
headeropt = {};
headeropt  = ft_setopt(headeropt, 'headerformat',   ft_getopt(cfg, 'headerformat'));        % is passed to low-level function, empty implies autodetection
headeropt  = ft_setopt(headeropt, 'readbids',       ft_getopt(cfg, 'readbids'));            % is passed to low-level function
headeropt  = ft_setopt(headeropt, 'coordsys',       ft_getopt(cfg, 'coordsys', 'head'));    % is passed to low-level function
headeropt  = ft_setopt(headeropt, 'coilaccuracy',   ft_getopt(cfg, 'coilaccuracy'));        % is passed to low-level function
headeropt  = ft_setopt(headeropt, 'coildeffile',    ft_getopt(cfg, 'coildeffile'));         % is passed to low-level function
headeropt  = ft_setopt(headeropt, 'checkmaxfilter', ft_getopt(cfg, 'checkmaxfilter'));      % this allows to read non-maxfiltered neuromag data recorded with internal active shielding
headeropt  = ft_setopt(headeropt, 'chantype',       ft_getopt(cfg, 'chantype', {}));        % 2017.10.10 AB required for NeuroOmega files
headeropt  = ft_setopt(headeropt, 'password',       ft_getopt(cfg, 'password'));            % this allows to read data from MED 1.0, MEF 3.0 and MEF 2.1 files
headeropt  = ft_setopt(headeropt, 'cache',          ft_getopt(cfg, 'cache'));

if ~isfield(cfg, 'feedback')
  if strcmp(cfg.method, 'channel')
    cfg.feedback = 'none';
  else
    cfg.feedback = 'text';
  end
end

% support for the following options was removed on 20 August 2004 in Revision 1.46
if isfield(cfg, 'emgchannel'), ft_error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'emghpfreq'),  ft_error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'emgrectify'), ft_error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'emghilbert'), ft_error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'eegchannel'), ft_error('EEG specific preprocessing is not supported any more'); end
if isfield(cfg, 'resamplefs'), ft_error('resampling is not supported any more, see RESAMPLEDATA'); end

if isfield(cfg, 'lnfilter') && strcmp(cfg.lnfilter, 'yes')
  ft_error('line noise filtering using the option cfg.lnfilter is not supported any more, use cfg.bsfilter instead')
end

% this relates to a previous fix to handle 32 bit neuroscan data
if isfield(cfg, 'nsdf')
  % FIXME this should be handled by ft_checkconfig, but ft_checkconfig does not allow yet for
  % specific errors in the case of forbidden fields
  ft_error('The use of cfg.nsdf is deprecated. FieldTrip tries to determine the bit resolution automatically. You can overrule this by specifying cfg.dataformat and cfg.headerformat. See: http://www.fieldtriptoolbox.org/faq/i_have_problems_reading_in_neuroscan_.cnt_files._how_can_i_fix_this');
end

if isfield(cfg, 'export') && ~isempty(cfg.export)
  % export the data to an output file
  if ~strcmp(cfg.method, 'trial')
    ft_error('exporting to an output file is only possible when processing all channels at once')
  end
end

% data can be passed by the user, but might also have been loaded from cfg.inputfile
hasdata = exist('data', 'var');

if hasdata
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do preprocessing of data that has already been read into memory
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % this is used to convert the data back to timelock later
  convert = ft_datatype(data);

  % check if the input data is valid for this function, the input data must be raw
  data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'hassampleinfo', 'yes');

  % check if the input cfg is valid for this function
  cfg = ft_checkconfig(cfg, 'forbidden',   {'trl', 'dataset', 'datafile', 'headerfile'});

  if cfg.padding>0
    if strcmp(cfg.dftfilter, 'yes') || ...
        strcmp(cfg.lpfilter, 'yes') || ...
        strcmp(cfg.hpfilter, 'yes') || ...
        strcmp(cfg.bpfilter, 'yes') || ...
        strcmp(cfg.bsfilter, 'yes') || ...
        strcmp(cfg.medianfilter, 'yes')
      padding = round(cfg.padding * data.fsample);
      if strcmp(cfg.padtype, 'data')
        ft_warning('datapadding not possible with in-memory data - padding will be performed by data mirroring');
        cfg.padtype = 'mirror';
      end
    else
      % no filtering will be done, hence no padding is necessary
      padding = 0;
    end
    % update the configuration (in seconds) for external reference
    cfg.padding = padding / data.fsample;
  else
    % no padding was requested
    padding = 0;
  end

  % some options don't make sense on component data
  if isfield(data, 'comp')
    if ~isempty(cfg.montage)
      ft_error('the application of a montage on component data is not supported');
    end
    if strcmp(cfg.reref, 'yes')
      ft_error('rereferencing component data is not supported');
    end
  end

  % set the defaults
  cfg.trials = ft_getopt(cfg, 'trials', 'all', 1);

  % select trials of interest
  tmpcfg = keepfields(cfg, {'trials', 'channel', 'latency', 'tolerance', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
  data   = ft_selectdata(tmpcfg, data);
  % restore the provenance information
  [cfg, data] = rollback_provenance(cfg, data);

  % this will contain the newly processed data
  % some fields from the input should be copied over in the output
  dataout = keepfields(data, {'hdr', 'fsample', 'grad', 'elec', 'opto', 'sampleinfo', 'trialinfo', 'topo', 'topolabel', 'unmixing'});

  ft_progress('init', cfg.feedback, 'preprocessing');
  ntrl = length(data.trial);
  dataout.trial = cell(1, ntrl);
  dataout.time  = cell(1, ntrl);
  for i=1:ntrl
    ft_progress(i/ntrl, 'preprocessing trial %d from %d\n', i, ntrl);
    nsamples = numel(data.time{i});

    % pad data by mirroring
    if nsamples>padding
      % the trial is already longer than the total length requested
      begpadding = 0;
      endpadding = 0;
      if padding > 0
        ft_warning('no padding applied because the padding duration is shorter than the trial');
      end
    else
      switch cfg.paddir
        case 'both'
          % begpadding+nsamples+endpadding = total length of data
          begpadding = ceil((padding-nsamples)/2);
          endpadding = floor((padding-nsamples)/2);
        case 'left'
          begpadding = padding-nsamples;
          endpadding = 0;
        case 'right'
          begpadding = 0;
          endpadding = padding-nsamples;
        otherwise
          ft_error('unsupported requested direction of padding');
      end
    end

    data.trial{i} = ft_preproc_padding(data.trial{i}, cfg.padtype, begpadding, endpadding);
    data.time{i}  = ft_preproc_padding(data.time{i}, 'nan',        begpadding, endpadding); % pad time-axis with nans (see bug2220)
    % do the filtering etc.
    [dataout.trial{i}, dataout.label, dataout.time{i}, cfg] = preproc(data.trial{i}, data.label,  data.time{i}, cfg, begpadding, endpadding);

  end % for all trials

  % convert back to input type if necessary
  switch convert
    case 'timelock'
      dataout = ft_checkdata(dataout, 'datatype', 'timelock');
    otherwise
      % keep the output as it is
  end
  ft_progress('close');

else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the data from file and do the preprocessing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if isfield(cfg, 'trialdef') && ~isfield(cfg, 'trl')
    ft_error('you must call FT_DEFINETRIAL prior to FT_PREPROCESSING');
  end

  % check if the input cfg is valid for this function
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
  cfg = ft_checkconfig(cfg, 'required',   {'headerfile', 'datafile'});
  cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
  cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

  % read the header
  hdr = ft_read_header(cfg.headerfile, headeropt{:});

  % this option relates to reading over trial boundaries in a pseudo-continuous dataset
  if ~isfield(cfg, 'continuous')
    if hdr.nTrials==1
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
  end

  if ~isfield(cfg, 'trl')
    % treat the data as continuous if possible, otherwise define all trials as indicated in the header
    if strcmp(cfg.continuous, 'yes')
      trl = zeros(1, 3);
      trl(1,1) = 1;
      trl(1,2) = hdr.nSamples*hdr.nTrials;
      trl(1,3) = -hdr.nSamplesPre;
    else
      trl = zeros(hdr.nTrials, 3);
      for i=1:hdr.nTrials
        trl(i,1) = (i-1)*hdr.nSamples + 1;
        trl(i,2) = (i  )*hdr.nSamples    ;
        trl(i,3) = -hdr.nSamplesPre;
      end
    end
    cfg.trl = trl;
  elseif ischar(cfg.trl)
    % load the trial information from file
    cfg.trl = loadvar(cfg.trl, 'trl');
  end

  % the code further down expects an Nx3 matrix with begsample, endsample and offset
  assert(size(cfg.trl,2)>=3, 'incorrect specification of cfg.trl');
  if istable(cfg.trl)
    trl = table2array(cfg.trl(:,1:3));
  else
    trl = cfg.trl(:,1:3);
  end

  % this should be a cell-array
  if ~iscell(cfg.channel) && ischar(cfg.channel)
    cfg.channel = {cfg.channel};
  end

  % this should be a cell-array
  if ~iscell(cfg.refchannel) && ischar(cfg.refchannel)
    cfg.refchannel = {cfg.refchannel};
  end

  % do a sanity check for the re-referencing
  if strcmp(cfg.reref, 'no') && ~isempty(cfg.refchannel)
    ft_warning('no re-referencing is performed');
    cfg.refchannel = {};
  end

  % translate the channel groups (like 'all' and 'MEG') into real labels
  cfg.channel = ft_channelselection(cfg.channel, hdr);
  assert(~isempty(cfg.channel), 'the selection of channels is empty');

  if ~isempty(cfg.implicitref)
    % add the label of the implicit reference channel to these cell-arrays
    cfg.channel = cat(1, cfg.channel(:), cfg.implicitref);
  end
  cfg.refchannel = ft_channelselection(cfg.refchannel, cfg.channel);

  % determine the length in samples to which the data should be padded before filtering is applied
  % the filter padding is done by reading a longer segment of data from the original data file
  if cfg.padding>0
    if strcmp(cfg.dftfilter, 'yes') || ...
        strcmp(cfg.lpfilter, 'yes') || ...
        strcmp(cfg.hpfilter, 'yes') || ...
        strcmp(cfg.bpfilter, 'yes') || ...
        strcmp(cfg.bsfilter, 'yes') || ...
        strcmp(cfg.medianfilter, 'yes')
      padding = round(cfg.padding * hdr.Fs);
    else
      % no filtering will be done, hence no padding is necessary
      padding = 0;
    end
    % update the configuration (in seconds) for external reference
    cfg.padding = padding / hdr.Fs;
  else
    % no padding was requested
    padding = 0;
  end

  if any(strmatch('reject',        fieldnames(cfg))) || ...
      any(strmatch('rejecteog',    fieldnames(cfg))) || ...
      any(strmatch('rejectmuscle', fieldnames(cfg))) || ...
      any(strmatch('rejectjump',   fieldnames(cfg)))
    % this is only for backward compatibility
    ft_error('you should call FT_REJECTARTIFACT prior to FT_PREPROCESSING, please update your scripts');
  end

  ntrl = size(trl,1);
  if ntrl<1
    ft_error('no trials were selected for preprocessing, see FT_DEFINETRIAL for help');
  end

  % compute the template for MCG and the QRS latency indices, and add it to the configuration
  if strcmp(cfg.removemcg, 'yes')
    cfg = template_mcg(cfg);
    mcgchannel = ft_channelselection(cfg.artfctdef.mcg.channel, hdr.label);
    mcgindx    = match_str(cfg.channel, mcgchannel);
    for i=1:length(mcgchannel)
      ft_info('removing mcg on channel %s\n', mcgchannel{i});
    end
  end

  % determine the channel numbers of interest for preprocessing
  [chnindx, rawindx] = match_str(cfg.channel, hdr.label);

  if strcmp(cfg.method, 'channel')
    % read one channel at a time, loop over channels and over trials
    chnloop = mat2cell(chnindx, ones(length(chnindx), 1), 1);
    rawloop = mat2cell(rawindx, ones(length(chnindx), 1), 1);

  elseif strcmp(cfg.method, 'trial')
    % read all channels simultaneously, only loop trials
    chnloop = {chnindx};
    rawloop = {rawindx};

  else
    ft_error('unsupported option for cfg.method');
  end

  for j=1:length(chnloop)
    % read one channel group at a time, this speeds up combined datasets
    % a multiplexed dataformat is faster if you read all channels, one trial at a time
    chnindx = chnloop{j};
    rawindx = rawloop{j};

    ft_info('processing channel { %s}\n', sprintf('''%s'' ', hdr.label{rawindx}));

    % initialize cell arrays
    cutdat = cell(1, ntrl);
    time   = cell(1, ntrl);

    ft_progress('init', cfg.feedback, 'reading and preprocessing');

    for i=1:ntrl
      ft_progress(i/ntrl, 'reading and preprocessing trial %d from %d\n', i, ntrl);
      % data padding is used for filtering and line noise removal
      nsamples = trl(i,2)-trl(i,1)+1;
      if nsamples>padding
        % the trial is already longer than the total length requested
        begsample  = trl(i,1);
        endsample  = trl(i,2);
        offset     = trl(i,3);
        begpadding = 0;
        endpadding = 0;
        if padding > 0
          ft_warning('no padding applied because the padding duration is shorter than the trial');
        end
      else
        switch cfg.paddir
          case 'both'
            % begpadding+nsamples+endpadding = total length of raw data that will be read
            begpadding = ceil((padding-nsamples)/2);
            endpadding = floor((padding-nsamples)/2);
          case 'left'
            begpadding = padding-nsamples;
            endpadding = 0;
          case 'right'
            begpadding = 0;
            endpadding = padding-nsamples;
          otherwise
            ft_error('unsupported requested direction of padding');
        end

        if strcmp(cfg.padtype, 'data')
          begsample  = trl(i,1) - begpadding;
          endsample  = trl(i,2) + endpadding;
        else
          % padding will be done below
          begsample  = trl(i,1);
          endsample  = trl(i,2);
        end
        if begsample<1
          ft_warning('cannot apply enough padding at begin of file');
          begpadding = begpadding - (1 - begsample);
          begsample  = 1;
        end
        if endsample>(hdr.nSamples*hdr.nTrials)
          ft_warning('cannot apply enough padding at end of file');
          endpadding = endpadding - (endsample - hdr.nSamples*hdr.nTrials);
          endsample  = hdr.nSamples*hdr.nTrials;
        end
        offset = trl(i,3) - begpadding;
      end

      % read the raw data with padding on both sides of the trial - this
      % includes datapadding
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', rawindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat, headeropt{:});

      % convert the data to another numeric precision, i.e. double, single or int32
      if ~isempty(cfg.precision)
        dat = cast(dat, cfg.precision);
      end

      % pad in case of no datapadding
      if ~strcmp(cfg.padtype, 'data')
        dat = ft_preproc_padding(dat, cfg.padtype, begpadding, endpadding);
        tim = offset2time(offset+begpadding, hdr.Fs, size(dat,2));
      else
        tim = offset2time(offset, hdr.Fs, size(dat,2));
      end

      % do the preprocessing on the padded trial data and remove the padding after filtering
      [cutdat{i}, label, time{i}, cfg] = preproc(dat, hdr.label(rawindx), tim, cfg, begpadding, endpadding);

      if isfield(cfg, 'export') && ~isempty(cfg.export)
        % write the processed data to an original manufacturer format file
        newhdr        = [];
        newhdr.Fs     = hdr.Fs;
        newhdr.label  = label;
        newhdr.nChans = length(newhdr.label);
        % only append for the second and consecutive trials
        ft_write_data(cfg.export.dataset, cutdat{i}, 'dataformat', cfg.export.dataformat, 'header', newhdr, 'append', i~=1);
        if nargout==0
          % don't keep the processed data in memory
          cutdat(i) = [];
        end
      end

    end % for all trials
    ft_progress('close');

    % don't keep hdr.orig in the output data if it is too large
    % hdr.orig can be large when caching data from specific file formats, such as bci2000_dat and mega_neurone
    if isfield(hdr, 'orig')
      s = hdr.orig;
      s = whos('s');
      if s.bytes>3*1024^2
        hdr = rmfield(hdr, 'orig');
      end
    end

    dataout                    = [];
    dataout.hdr                = hdr;                 % header details of the datafile
    dataout.label              = label;               % labels of channels that have been read, can be different from labels in file due to montage
    dataout.time               = time;                % vector with the timeaxis for each individual trial
    dataout.trial              = cutdat;
    dataout.fsample            = hdr.Fs;
    if istable(cfg.trl)
      % we always want the sampleinfo to be numeric
      dataout.sampleinfo       = table2array(cfg.trl(:,1:2));
    else
      dataout.sampleinfo       = cfg.trl(:,1:2);
    end
    if size(cfg.trl,2) > 3
      dataout.trialinfo        = cfg.trl(:,4:end);    % this can be a numeric array or a table
    end
    if isfield(hdr, 'grad')
      dataout.grad             = hdr.grad;            % MEG gradiometer information in header (f.e. headerformat = 'ctf_ds')
    end
    if isfield(hdr, 'elec')
      dataout.elec             = hdr.elec;            % EEG electrode information in header (f.e. headerformat = 'neuromag_fif')
    end
    if isfield(hdr, 'opto')
      dataout.opto             = hdr.opto;            % NIRS optode information in header (f.e. headerformat = 'artinis')
    end

  end % for all channel groups

end % if hasdata

if strcmp(cfg.updatesens, 'yes')
  % updating the sensor descriptions can be done on basis of the montage or the rereference settings
  if ~isempty(cfg.montage) && ~isequal(cfg.montage, 'no')
    montage = cfg.montage;
  elseif strcmp(cfg.reref, 'yes')
    if strcmp(cfg.refmethod, 'bipolar') || strcmp(cfg.refmethod, 'avg') || strcmp(cfg.refmethod, 'laplace')
      tmpcfg = keepfields(cfg, {'refmethod', 'implicitref', 'refchannel', 'channel', 'groupchans'});
      tmpcfg.showcallinfo = 'no';
      montage = ft_prepare_montage(tmpcfg, data);
    else
      % do not update the sensor description
      montage = [];
    end
  else
    % do not update the sensor description
    montage = [];
  end

  if ~isempty(montage)
    % apply the linear projection also to the sensor description
    if issubfield(montage, 'type')
      bname = montage.type;
    else
      bname = 'preproc';
    end
    if isfield(dataout, 'grad')
      ft_info('applying the montage to the grad structure\n');
      dataout.grad = ft_apply_montage(dataout.grad, montage, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
    end
    if isfield(dataout, 'elec')
      ft_info('applying the montage to the elec structure\n');
      dataout.elec = ft_apply_montage(dataout.elec, montage, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
    end
    if isfield(dataout, 'opto')
      ft_info('applying the montage to the opto structure\n');
      dataout.opto = ft_apply_montage(dataout.opto, montage, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
    end
  end
end % if updatesens

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous data

% rename the output variable to accommodate the savevar postamble
data = dataout;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
