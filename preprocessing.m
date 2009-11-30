function [data] = preprocessing(cfg, data);

% PREPROCESSING reads MEG and/or EEG data according to user-specified trials
% and applies several user-specified preprocessing steps to the signals.
%
% Use as
%   [data] = preprocessing(cfg)
% or
%   [data] = preprocessing(cfg, data)
%
% The first input argument "cfg" is the configuration structure, which
% contains all details for the dataset filenames, trials and the
% preprocessing options. You can only do preprocessing after defining the
% segments of data to be read from the file (i.e. the trials), which is for
% example done based on the occurence of a trigger in the data.
%
% If you are calling PREPROCESSING with only the configuration as first
% input argument and the data still has to be read from file, you should
% specify
%   cfg.dataset      = string with the filename
%   cfg.trl          = Nx3 matrix with the trial definition, see DEFINETRIAL
%   cfg.padding      = length to which the trials are padded for filtering
%   cfg.continuous   = 'yes' or 'no' whether the file contains continuous data
%                      (default is determined automatic)
%
% Instead of specifying the dataset, you can also explicitely specify the
% name of the file containing the header information and the name of the
% file containing the data, using
%   cfg.datafile     = string with the filename
%   cfg.headerfile   = string with the filename
%
% If you are calling PREPROCESSING with also the second input argument
% "data", then that should contain data that was already read from file in
% a previous call to PREPROCESSING. In that case only the configuration
% options below apply.
%
% The channels that will be read and/or preprocessed are specified with
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see CHANNELSELECTION for details
%
% The preprocessing options for the selected channels are specified with
%   cfg.lpfilter      = 'no' or 'yes'  lowpass filter
%   cfg.hpfilter      = 'no' or 'yes'  highpass filter
%   cfg.bpfilter      = 'no' or 'yes'  bandpass filter
%   cfg.bsfilter      = 'no' or 'yes'  bandstop filter
%   cfg.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   cfg.medianfilter  = 'no' or 'yes'  jump preserving median filter
%   cfg.lpfreq        = lowpass  frequency in Hz
%   cfg.hpfreq        = highpass frequency in Hz
%   cfg.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.bsfreq        = bandstop frequency range, specified as [low high] in Hz
%   cfg.dftfreq       = line noise frequencies for DFT filter, default [50 100 150] Hz
%   cfg.lpfiltord     = lowpass  filter order
%   cfg.hpfiltord     = highpass filter order
%   cfg.bpfiltord     = bandpass filter order
%   cfg.bsfiltord     = bandstop filter order
%   cfg.lpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.hpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.bpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.bsfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.lpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.hpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.bpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.bsfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.medianfiltord = length of median filter
%   cfg.blc           = 'no' or 'yes'
%   cfg.blcwindow     = [begin end] in seconds, the default is the complete trial
%   cfg.detrend       = 'no' or 'yes', this is done on the complete trial
%   cfg.polyremoval   = 'no' or 'yes', this is done on the complete trial
%   cfg.polyorder     = polynome order (default = 2)
%   cfg.derivative    = 'no' (default) or 'yes', computes the first order derivative of the data
%   cfg.hilbert       = 'no', 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag' or 'angle' (default = 'no')
%   cfg.rectify       = 'no' or 'yes'
%   cfg.precision     = 'single' or 'double' (default = 'double')
%
% Preprocessing options that you should only use for EEG data are
%   cfg.reref         = 'no' or 'yes' (default = 'no')
%   cfg.refchannel    = cell-array with new EEG reference channel(s)
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%   cfg.montage       = 'no' or a montage structure (default = 'no')
%
% Preprocessing options that you should only use when you are calling PREPROCESSING with
% also the second input argument "data" are
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')

% Undocumented local options:
% cfg.artfctdef
% cfg.removemcg
% cfg.output.dataset
% cfg.output.dataformat
%
% This function depends on PREPROC which has the following options:
% cfg.absdiff
% cfg.boxcar
% cfg.polyremoval, documented
% cfg.polyorder, documented
% cfg.blc, documented
% cfg.blcwindow, documented
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

% Copyright (C) 2003-2007, Robert Oostenveld, SMI, FCDC
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'channel'),      cfg.channel = {'all'};        end
if ~isfield(cfg, 'removemcg'),    cfg.removemcg = 'no';         end
if ~isfield(cfg, 'removeeog'),    cfg.removeeog = 'no';         end
if ~isfield(cfg, 'feedback'),     cfg.feedback = 'text';        end
if ~isfield(cfg, 'precision'),    cfg.precision = 'double';     end
if ~isfield(cfg, 'padding'),      cfg.padding = 0;              end % padding is only done when filtering
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

if nargin>1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do preprocessing of data that has already been read
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % the input data must be raw
  data = checkdata(data, 'datatype', 'raw', 'hasoffset', 'yes');
  
  % check if the input cfg is valid for this function
  cfg = checkconfig(cfg, 'forbidden',   {'trl', 'dataset', 'datafile', 'headerfile'});
  
  if cfg.padding>0
    error('cfg.padding should be zero, since filter padding is only possible while reading the data from file');
  end
  
  % set the defaults
  if ~isfield(cfg, 'trials'), cfg.trials = 'all'; end
  
  % select trials of interest
  if ~strcmp(cfg.trials, 'all')
    if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
    fprintf('selecting %d trials\n', length(cfg.trials));
    data.trial  = data.trial(cfg.trials);
    data.time   = data.time(cfg.trials);
    data.offset = data.offset(cfg.trials);
  end
  
  % translate the channel groups (like 'all' and 'MEG') into real labels
  cfg.channel = channelselection(cfg.channel, data.label);
  rawindx = match_str(data.label, cfg.channel);
  
  % this will contain the newly processed data
  dataout = [];
  
  progress('init', cfg.feedback, 'preprocessing');
  ntrl = length(data.trial);
  for i=1:ntrl
    progress(i/ntrl, 'preprocessing trial %d from %d\n', i, ntrl);
    % do the preprocessing on the selected channels
    [dataout.trial{i}, dataout.label, dataout.time{i}, cfg] = preproc(data.trial{i}(rawindx,:), data.label(rawindx), data.fsample, cfg, data.offset(i));
    
    if isfield(cfg, 'output') && ~isempty(cfg.output)
      % write the processed data to file
      newhdr        = [];
      newhdr.Fs     = hdr.Fs;
      newhdr.label  = label;
      newhdr.nChans = length(newhdr.label);
      
      % only append for the second and consecutive trials
      write_data(cfg.output.dataset, dataout.trial{i}, 'dataformat', cfg.output.dataformat, 'header', newhdr, 'append', i~=1);
      
      if nargout==0
        % don't keep the data in memory
        dataout.trial{i} = [];
      end
    end
    
  end % for all trials
  progress('close');
  
  % update the trial definition (trl) in case of trial selection
  if ~strcmp(cfg.trials, 'all')
    % try to locate the trl in the nested configuration
    if isfield(data, 'cfg')
      trl = findcfg(data.cfg, 'trl');
    else
      trl = [];
    end
    if isempty(trl)
      % a trial definition is expected in each continuous data set
      warning('could not locate the trial definition ''trl'' in the data structure');
    else
      cfg.trlold=trl;
      cfg.trl=trl(cfg.trials,:);
    end
  end
  
  % get the output cfg
  cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');
  
  % remember the configuration details of the input data
  if isfield(data, 'cfg'); cfg.previous = data.cfg; end
  
  % take along relevant fields of input data to output data
  if isfield(data, 'hdr'),      dataout.hdr     = data.hdr;         end
  if isfield(data, 'fsample'),  dataout.fsample = data.fsample;     end
  if isfield(data, 'grad'),     dataout.grad    = data.grad;        end
  if isfield(data, 'elec'),     dataout.elec    = data.elec;        end
  if isfield(data, 'offset'),   dataout.offset  = data.offset;      end

  % replace the input data with the output data
  data = dataout;
  
else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the data from file and do the preprocessing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if isfield(cfg, 'trialdef') && ~isfield(cfg, 'trl')
    error('you must call DEFINETRIAL prior to PREPROCESSING');
  end
  
  % check if the input cfg is valid for this function
  cfg = checkconfig(cfg, 'dataset2files', {'yes'});
  cfg = checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  cfg = checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
  cfg = checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});
  
  % read the header
  hdr = read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
  
  % this option relates to reading over trial boundaries in a pseudo-continuous dataset
  if ~isfield(cfg, 'continuous')
    if hdr.nTrials==1
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
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
  cfg.channel = channelselection(cfg.channel, hdr.label);
  
  if ~isempty(cfg.implicitref)
    % add the label of the implicit reference channel to these cell-arrays
    cfg.channel    = cat(1, cfg.channel(:), cfg.implicitref);
  end
  cfg.refchannel = channelselection(cfg.refchannel, cfg.channel);
  
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
  
  if any(strmatch('reject',       fieldnames(cfg))) || ...
      any(strmatch('rejecteog',    fieldnames(cfg))) || ...
      any(strmatch('rejectmuscle', fieldnames(cfg))) || ...
      any(strmatch('rejectjump',   fieldnames(cfg)))
    % this is only for backward compatibility
    error('you should call REJECTARTIFACT prior to PREPROCESSING, please update your scripts');
  end
  
  ntrl = size(cfg.trl,1);
  if ntrl<1
    error('no trials were selected for preprocessing, see DEFINETRIAL for help');
  end
  
  % compute the template for MCG and the QRS latency indices, and add it to the configuration
  if strcmp(cfg.removemcg, 'yes')
    cfg = template_mcg(cfg);
    mcgchannel = channelselection(cfg.artfctdef.mcg.channel, hdr.label);
    mcgindx    = match_str(cfg.channel, mcgchannel);
    for i=1:length(mcgchannel)
      fprintf('removing mcg on channel %s\n', mcgchannel{i});
    end
  end
  
  % determine the channel numbers of interest for preprocessing
  [chnindx, rawindx] = match_str(cfg.channel, hdr.label);
  
  progress('init', cfg.feedback, 'reading and preprocessing');
  for i=1:ntrl
    progress(i/ntrl, 'reading and preprocessing trial %d from %d\n', i, ntrl);
    % non-zero padding is used for filtering and line noise removal
    nsamples = cfg.trl(i,2)-cfg.trl(i,1)+1;
    if nsamples>padding
      % the trial is already longer than the total lenght requested
      begsample  = cfg.trl(i,1);
      endsample  = cfg.trl(i,2);
      begpadding = 0;
      endpadding = 0;
    else
      % begpadding+nsamples+endpadding = total length of raw data that will be read
      begpadding = ceil((padding-nsamples)/2);
      endpadding = floor((padding-nsamples)/2);
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
    end
    
    % ONLY RELEVANT FOR NEUROSCAN CNT
    if ~isfield(cfg, 'nsdf')
      hdr.nsdf=16;
    else
      hdr.nsdf=cfg.nsdf;
    end
    
    % read the raw data with padding on both sides of the trial
    dat = read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', rawindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
    
    % do the preprocessing on the padded trial data and remove the padding after filtering
    [cutdat{i}, label, time{i}, cfg] = preproc(dat, hdr.label(rawindx), hdr.Fs, cfg, cfg.trl(i,3), begpadding, endpadding);
    
    if isfield(cfg, 'output') && ~isempty(cfg.output)
      % write the processed data to file
      newhdr        = [];
      newhdr.Fs     = hdr.Fs;
      newhdr.label  = label;
      newhdr.nChans = length(newhdr.label);
      
      % only append for the second and consecutive trials
      write_data(cfg.output.dataset, cutdat{i}, 'dataformat', cfg.output.dataformat, 'header', newhdr, 'append', i~=1);
      
      if nargout==0
        % don't keep the data in memory
        cutdat{i} = [];
      end
    end
    
    % ONLY RELEVANT FOR NEUROSCAN CNT
    hdr=rmfield(hdr,'nsdf');
    
  end % for all trials
  progress('close');
  
  % collect the results
  data.hdr                = hdr;                  % header details of the datafile
  data.label              = label;                % labels of channels that have been read, can be different from labels in file due to montage
  data.trial              = cutdat;               % cell-array with TIMExCHAN
  data.time               = time;                 % vector with the timeaxis for each individual trial
  data.fsample            = hdr.Fs;
  if isfield(hdr, 'grad')
    data.grad             = hdr.grad;             % gradiometer system in head coordinates
  end
  
  % get the output cfg
  cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');
  
end % if nargin>1

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id$';
% remember the exact configuration details in the output
data.cfg = cfg;

