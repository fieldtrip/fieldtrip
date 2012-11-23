function [data] = ft_preprocessing(cfg, data)

% FT_PREPROCESSING reads MEG and/or EEG data according to user-specified trials
% and applies several user-specified preprocessing steps to the signals.
%
% Use as
%   [data] = ft_preprocessing(cfg)
% or
%   [data] = ft_preprocessing(cfg, data)
%
% The first input argument "cfg" is the configuration structure, which
% contains all details for the dataset filenames, trials and the
% preprocessing options. You can only do preprocessing after defining the
% segments of data to be read from the file (i.e. the trials), which is for
% example done based on the occurence of a trigger in the data.
%
% If you are calling FT_PREPROCESSING with only the configuration as first
% input argument and the data still has to be read from file, you should
% specify
%   cfg.dataset      = string with the filename
%   cfg.trl          = Nx3 matrix with the trial definition, see FT_DEFINETRIAL
%   cfg.padding      = length to which the trials are padded for filtering (default = 0)
%   cfg.padtype      = string, type of padding (default: 'data' padding or
%                      'mirror', depending on feasibility)
%   cfg.continuous   = 'yes' or 'no' whether the file contains continuous data
%                      (default is determined automatic)
%
% Instead of specifying the dataset, you can also explicitely specify the
% name of the file containing the header information and the name of the
% file containing the data, using
%   cfg.datafile     = string with the filename
%   cfg.headerfile   = string with the filename
%
% If you are calling FT_PREPROCESSING with also the second input argument
% "data", then that should contain data that was already read from file in
% a previous call to FT_PREPROCESSING. In that case only the configuration
% options below apply.
%
% The channels that will be read and/or preprocessed are specified with
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
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
%   cfg.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.bsfreq        = bandstop frequency range, specified as [low high] in Hz
%   cfg.dftfreq       = line noise frequencies in Hz for DFT filter (default = [50 100 150])
%   cfg.lpfiltord     = lowpass  filter order (default set in low-level function)
%   cfg.hpfiltord     = highpass filter order (default set in low-level function)
%   cfg.bpfiltord     = bandpass filter order (default set in low-level function)
%   cfg.bsfiltord     = bandstop filter order (default set in low-level function)
%   cfg.lpfilttype    = digital filter type, 'but' or 'fir' or 'firls' (default = 'but')
%   cfg.hpfilttype    = digital filter type, 'but' or 'fir' or 'firls' (default = 'but')
%   cfg.bpfilttype    = digital filter type, 'but' or 'fir' or 'firls' (default = 'but')
%   cfg.bsfilttype    = digital filter type, 'but' or 'fir' or 'firls' (default = 'but')
%   cfg.lpfiltdir     = filter direction, 'twopass', 'onepass' or 'onepass-reverse' (default = 'twopass') 
%   cfg.hpfiltdir     = filter direction, 'twopass', 'onepass' or 'onepass-reverse' (default = 'twopass') 
%   cfg.bpfiltdir     = filter direction, 'twopass', 'onepass' or 'onepass-reverse' (default = 'twopass') 
%   cfg.bsfiltdir     = filter direction, 'twopass', 'onepass' or 'onepass-reverse' (default = 'twopass') 
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
%
% Preprocessing options that you should only use for EEG data are
%   cfg.reref         = 'no' or 'yes' (default = 'no')
%   cfg.refchannel    = cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%   cfg.montage       = 'no' or a montage structure (default = 'no')
%
% Preprocessing options that you should only use when you are calling FT_PREPROCESSING with
% also the second input argument "data" are
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%
% Preprocessing options that you should only use when you are calling
% FT_PREPROCESSING with a single cfg input argument are
%   cfg.method        = 'trial' or 'channel', read data per trial or per channel (default = 'trial')
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_DEFINETRIAL, FT_REDEFINETRIAL, FT_APPENDDATA, FT_APPENDSPIKE

% Guidelines for use in an analysis pipeline: 
% After FT_PREPROCESSING you will have raw data represented as a single 
% continuous segment or as multiple data segments that often correspond to 
% trials in an experiment.
% This usually serves as input for one of the following functions:
%    * FT_TIMELOCKANALYSIS  to compute event-related fields or potentials
%    * FT_FREQANALYSIS      to compute the frequency or time-frequency representation
%    * FT_PREPROCESSING     if you want to apply additional temporal filters, baseline correct, rereference or apply an EEG montage
%    * FT_APPENDDATA        if you have preprocessed separate conditions or datasets and want to combine them
%    * FT_REDEFINETRIAL     if you want to cut the data segments into smaller pieces or want to change the time axes
%    * FT_DATABROWSER       to inspect the data and check for artifacts
%    * FT_REJECTVISUAL      to inspect the data and remove trials that contain artifacts
%    * FT_COMPONENTANALYSIS if you want to use ICA to remove artifacts

% Undocumented local options:
% cfg.paddir = direction of padding, 'left'/'right'/'both' (default = 'both')
% cfg.artfctdef
% cfg.removemcg
% cfg.montage (in combination with meg-data in the input) applies montage
%              to both data and grad-structure)
% You can use this function to read data from one format, filter it, and
% write it to disk in another format. The reading is done either as one
% long continuous segment or in multiple trials. This is achieved by
%   cfg.export.dataset    = string with the output file name
%   cfg.export.dataformat = string describing the output file format, see FT_WRITE_DATA

% This function depends on PREPROC which has the following options:
% cfg.absdiff
% cfg.boxcar
% cfg.polyremoval, documented
% cfg.polyorder, documented
% cfg.demean, documented
% cfg.baselinewindow, documented
% cfg.bpfilter, documented
% cfg.bpfiltord, documented
% cfg.bpfilttype, documented
% cfg.bpfreq, documented
% cfg.bsfilter, documented
% cfg.bsfiltord, documented
% cfg.bsfilttype, documented
% cfg.bsfreq, documented
% cfg.derivative, documented
% cfg.detrend, documented
% cfg.dftfilter, documented
% cfg.dftfreq, documented
% cfg.hilbert, documented
% cfg.hpfilter, documented
% cfg.hpfiltord, documented
% cfg.hpfilttype, documented
% cfg.hpfreq, documented
% cfg.implicitref, documented
% cfg.lpfilter, documented
% cfg.lpfiltord, documented
% cfg.lpfilttype, documented
% cfg.lpfreq, documented
% cfg.medianfilter, documented
% cfg.medianfiltord, documented
% cfg.rectify, documented
% cfg.refchannel, documented
% cfg.reref, documented

% Copyright (C) 2003-2012, Robert Oostenveld, SMI, FCDC
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig
ft_preamble loadvar data

% return immediately after distributed execution
if ~isempty(ft_getopt(cfg, 'distribute'))
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed', {'blcwindow', 'baselinewindow'});
cfg = ft_checkconfig(cfg, 'renamed', {'output', 'export'});

% set the defaults
if ~isfield(cfg, 'method'),       cfg.method = 'trial';         end
if ~isfield(cfg, 'channel'),      cfg.channel = {'all'};        end
if ~isfield(cfg, 'removemcg'),    cfg.removemcg = 'no';         end
if ~isfield(cfg, 'removeeog'),    cfg.removeeog = 'no';         end

if ~isfield(cfg, 'feedback'),
  if strcmp(cfg.method, 'channel')
    cfg.feedback = 'none';
  else
    cfg.feedback = 'text';
  end
end

if ~isfield(cfg, 'precision'),    cfg.precision = 'double';     end
if ~isfield(cfg, 'padding'),      cfg.padding = 0;              end % padding is only done when filtering
if ~isfield(cfg, 'paddir'),       cfg.paddir  = 'both';         end
if ~isfield(cfg, 'headerformat'), cfg.headerformat = [];        end % is passed to low-level function, empty implies autodetection
if ~isfield(cfg, 'dataformat'),   cfg.dataformat = [];          end % is passed to low-level function, empty implies autodetection

% these options relate to the actual preprocessing, it is neccessary to specify here because of padding
if ~isfield(cfg, 'dftfilter'),    cfg.dftfilter = 'no';         end
if ~isfield(cfg, 'lpfilter'),     cfg.lpfilter = 'no';          end
if ~isfield(cfg, 'hpfilter'),     cfg.hpfilter = 'no';          end
if ~isfield(cfg, 'bpfilter'),     cfg.bpfilter = 'no';          end
if ~isfield(cfg, 'bsfilter'),     cfg.bsfilter = 'no';          end
if ~isfield(cfg, 'medianfilter'), cfg.medianfilter = 'no';      end
% these options relate to the actual preprocessing, it is neccessary to specify here because of channel selection
if ~isfield(cfg, 'reref'),        cfg.reref = 'no';             end
if ~isfield(cfg, 'refchannel'),   cfg.refchannel = {};          end
if ~isfield(cfg, 'implicitref'),  cfg.implicitref = [];         end

cfg.padtype = ft_getopt(cfg, 'padtype', 'data');

% support for the following options was removed on 20 August 2004 in Revision 1.46
if isfield(cfg, 'emgchannel'), error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'emghpfreq'),  error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'emgrectify'), error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'emghilbert'), error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'eegchannel'), error('EEG specific preprocessing is not supported any more'); end
if isfield(cfg, 'resamplefs'), error('resampling is not supported any more, see RESAMPLEDATA'); end

if isfield(cfg, 'lnfilter') && strcmp(cfg.lnfilter, 'yes')
  error('line noise filtering using the option cfg.lnfilter is not supported any more, use cfg.bsfilter instead')
end

% this relates to a previous fix to handle 32 bit neuroscan data
if isfield(cfg, 'nsdf'),
  % FIXME this should be handled by ft_checkconfig, but ft_checkconfig does not allow yet for
  % specific errors in the case of forbidden fields
  error('the use of cfg.nsdf is deprecated. fieldtrip tries to determine the bit resolution automatically. you can overrule this by specifying cfg.dataformat and cfg.headerformat. see: http://fieldtrip.fcdonders.nl/faq/i_have_problems_reading_in_neuroscan_.cnt_files._how_can_i_fix_this');
end

if isfield(cfg, 'export') && ~isempty(cfg.export)
  % export the data to an output file
  if ~strcmp(cfg.method, 'trial')
    error('exporting to an output file is only possible when processing all channels at once')
  end
end

hasdata = exist('data', 'var');
if hasdata
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do preprocessing of data that has already been read into memory
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  convert = ft_datatype(data);
  % the input data must be raw
  data = ft_checkdata(data, 'datatype', 'raw', 'hassampleinfo', 'yes');

  % check if the input cfg is valid for this function
  cfg = ft_checkconfig(cfg, 'forbidden',   {'trl', 'dataset', 'datafile', 'headerfile'});

  if cfg.padding>0
    if strcmp(cfg.dftfilter, 'yes') || ...
        strcmp(cfg.lpfilter, 'yes') || ...
        strcmp(cfg.hpfilter, 'yes') || ...
        strcmp(cfg.bpfilter, 'yes') || ...
        strcmp(cfg.bsfilter, 'yes') || ...
        strcmp(cfg.medianfilter, 'yes')
      padding = round(cfg.padding * data.Fs);
      if strcmp(cfg.padtype, 'data')
        warning_once('datapadding not possible with in-memory data - padding will be performed by data mirroring');
        cfg.padtype = 'mirror';
      end
    else
      % no filtering will be done, hence no padding is neccessary
      padding = 0;
    end
    % update the configuration (in seconds) for external reference
    cfg.padding = padding / data.Fs;
  else
    % no padding was requested
    padding = 0;
  end  

  % set the defaults
  if ~isfield(cfg, 'trials'), cfg.trials = 'all'; end

  % select trials of interest
  if ~strcmp(cfg.trials, 'all')
    data = ft_selectdata(data, 'rpt', cfg.trials);
  end

  % translate the channel groups (like 'all' and 'MEG') into real labels
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  rawindx = match_str(data.label, cfg.channel);

  % this will contain the newly processed data
  dataout = [];
  % take along relevant fields of input data to output data
  if isfield(data, 'hdr'),      dataout.hdr     = data.hdr;         end
  if isfield(data, 'fsample'),  dataout.fsample = data.fsample;     end
  if isfield(data, 'grad'),     dataout.grad    = data.grad;        end
  if isfield(data, 'elec'),     dataout.elec    = data.elec;        end
  if isfield(data, 'sampleinfo'),  dataout.sampleinfo  = data.sampleinfo;  end
  if isfield(data, 'trialinfo'), dataout.trialinfo = data.trialinfo; end
  
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
          error('unsupported requested direction of padding');
      end
    end
          
    data.trial{i} = ft_preproc_padding(data.trial{i}, padtype, begpadding, endpadding);
        
    % do the preprocessing on the selected channels
    [dataout.trial{i}, dataout.label, dataout.time{i}, cfg] = preproc(data.trial{i}(rawindx,:), data.label(rawindx), data.time{i}, cfg, begpadding, endpadding);
    
  end % for all trials
  
  if isfield(dataout, 'grad') && isfield(cfg, 'montage') && ~strcmp(cfg.montage, 'no') && isstruct(cfg.montage)
    % apply the montage also to the MEG-sensor description
    if isfield(cfg.montage, 'type'),
      bname = cfg.montage.type;
    else
      bname = 'preproc';
    end
    dataout.grad = ft_apply_montage(dataout.grad, cfg.montage, 'feedback', 'none', 'keepunused', 'yes', 'balancename', bname);
  end
  
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
    error('you must call FT_DEFINETRIAL prior to FT_PREPROCESSING');
  end

  % check if the input cfg is valid for this function
  cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});
  cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
  cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

  % read the header
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);

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
  end

  % this should be a cell array
  if ~iscell(cfg.channel) && ischar(cfg.channel)
    cfg.channel = {cfg.channel};
  end

  % this should be a cell array
  if ~iscell(cfg.refchannel) && ischar(cfg.refchannel)
    cfg.refchannel = {cfg.refchannel};
  end

  % do a sanity check for the re-referencing
  if strcmp(cfg.reref, 'no') && ~isempty(cfg.refchannel)
    warning('no re-referencing is performed');
    cfg.refchannel = {};
  end

  % translate the channel groups (like 'all' and 'MEG') into real labels
  cfg.channel = ft_channelselection(cfg.channel, hdr.label);

  if ~isempty(cfg.implicitref)
    % add the label of the implicit reference channel to these cell-arrays
    cfg.channel    = cat(1, cfg.channel(:), cfg.implicitref);
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
      % no filtering will be done, hence no padding is neccessary
      padding = 0;
    end
    % update the configuration (in seconds) for external reference
    cfg.padding = padding / hdr.Fs;
  else
    % no padding was requested
    padding = 0;
  end

  if any(strmatch('reject',       fieldnames(cfg))) || ...
      any(strmatch('rejecteog',    fieldnames(cfg))) || ...
      any(strmatch('rejectmuscle', fieldnames(cfg))) || ...
      any(strmatch('rejectjump',   fieldnames(cfg)))
    % this is only for backward compatibility
    error('you should call FT_REJECTARTIFACT prior to FT_PREPROCESSING, please update your scripts');
  end

  ntrl = size(cfg.trl,1);
  if ntrl<1
    error('no trials were selected for preprocessing, see FT_DEFINETRIAL for help');
  end

  % compute the template for MCG and the QRS latency indices, and add it to the configuration
  if strcmp(cfg.removemcg, 'yes')
    cfg = template_mcg(cfg);
    mcgchannel = ft_channelselection(cfg.artfctdef.mcg.channel, hdr.label);
    mcgindx    = match_str(cfg.channel, mcgchannel);
    for i=1:length(mcgchannel)
      fprintf('removing mcg on channel %s\n', mcgchannel{i});
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
    error('unsupported option for cfg.method');
  end

  for j=1:length(chnloop)
    % read one channel group at a time, this speeds up combined datasets
    % a multiplexed dataformat is faster if you read all channels, one trial at a time
    chnindx = chnloop{j};
    rawindx = rawloop{j};

    fprintf('processing channel { %s}\n', sprintf('''%s'' ', hdr.label{rawindx}));
    
    % initialize cell arrays
    cutdat = cell(1, ntrl);
    time   = cell(1, ntrl);
    
    ft_progress('init', cfg.feedback, 'reading and preprocessing');
    
    for i=1:ntrl
      ft_progress(i/ntrl, 'reading and preprocessing trial %d from %d\n', i, ntrl);
      % non-zero padding is used for filtering and line noise removal
      nsamples = cfg.trl(i,2)-cfg.trl(i,1)+1;
      if nsamples>padding
        % the trial is already longer than the total length requested
        begsample  = cfg.trl(i,1);
        endsample  = cfg.trl(i,2);
        offset     = cfg.trl(i,3);
        begpadding = 0;
        endpadding = 0;
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
          error('unsupported requested direction of padding');
        end
        
        if strcmp(cfg.padding, 'data');
          begsample  = cfg.trl(i,1) - begpadding;
          endsample  = cfg.trl(i,2) + endpadding;
          if begsample<1
            warning('cannot apply enough padding at begin of file');
            begpadding = begpadding - (1 - begsample);
            begsample  = 1;
          end
          if endsample>(hdr.nSamples*hdr.nTrials)
            warning('cannot apply enough padding at end of file');
            endpadding = endpadding - (endsample - hdr.nSamples*hdr.nTrials);
            endsample  = hdr.nSamples*hdr.nTrials;
          end
          offset = cfg.trl(i,3) - begpadding;
        end
      end

      % read the raw data with padding on both sides of the trial
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', rawindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
      tim = offset2time(offset, hdr.Fs, size(dat,2));
      
      if ~strcmp(cfg.padtype, 'data')
        dat = ft_preproc_padding(dat, padtype, begpadding, endpadding);
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
          % don't keep th eprocessed data in memory
          cutdat(i) = [];
        end
      end

    end % for all trials
    ft_progress('close');

    dataout                    = [];
    dataout.hdr                = hdr;                  % header details of the datafile
    dataout.label              = label;                % labels of channels that have been read, can be different from labels in file due to montage
    dataout.time               = time;                 % vector with the timeaxis for each individual trial
    dataout.trial              = cutdat;
    dataout.fsample            = hdr.Fs;
    dataout.sampleinfo         = cfg.trl(:,1:2);
    if size(cfg.trl,2) > 3
        dataout.trialinfo      = cfg.trl(:,4:end);
    end
    if isfield(hdr, 'grad')
      dataout.grad             = hdr.grad;             % gradiometer system in head coordinates
    end

  end % for all channel groups

end % if hasdata

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble provenance

if hasdata
  ft_postamble previous data
end

% rename the output variable to accomodate the savevar postamble
data = dataout;

ft_postamble history data
ft_postamble savevar data

