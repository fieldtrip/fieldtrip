function [cfg] = spikedownsample(cfg)

% SPIKEDOWNSAMPLE takes electrophysiological data that was continuoudly
% sampled at 32KHz and preprocesses and downsamples it to obtain the LFP
% data, which can subsequently be processed in more detail.
%
% Use as
%   [cfg] = spikedownsample(cfg)
%
% The configuration should contain
%   cfg.dataset             = string with the input dataset
%   cfg.output              = string with the output dataset (default is determined automatic)
%   cfg.dataformat          = string with the output dataset format, see WRITE_DATA
%   cfg.channel             = Nx1 cell-array with selection of channels (default = 'all'),
%                             see CHANNELSELECTION for details
%   cfg.fsample             = desired sampling frequency in Hz (default = 1000)
%   cfg.method              = resampling method, can be 'resample', 'decimate' or 'downsample'
%   cfg.timestampdefinition = 'orig' or 'sample'
%   cfg.channelprefix       = string, will be added to channel name, e.g. 'lfp' -> 'lfp_ncs001' (default = [])
%   cfg.calibration         = optional scaling factor to apply to the data to convert it in uV, see below
%
% The Neuralynx acquisition system at the FCDC in Nijmegen makes use of
% Plexon headstages which have a amplification of 20x. The data that is
% written by the Neuralynx acquisition software therefore is 20x larger
% than the true microvolt values. When operating SPIKEDOWNSAMPLE on the
% *.ncs files that are recorded with the Neuralynx Cheetah software, the
% calibration should be set to 1/20. The raw dma file (*.nrd) and the
% splitted DMA files contains AD values that are not scaled in uV and
% require an additional factor of 64x. If you operate SPIKEDOWNSAMPLE  on
% raw dma files or on splitted DMA files, the calibration should be set to
% 1/(64*20).
%
% The default is to process the full dataset. You can select a latency range with
%   cfg.latency          = [begin end], default is [0 inf]
% or you can specify multiple latency segments with
%   cfg.latency          = [b1 e1; b2 e2; ...]
%
% Furthermore, the configuration can contain the following preprocessing options
%   cfg.preproc.lpfilter      = 'no' or 'yes'  lowpass filter
%   cfg.preproc.hpfilter      = 'no' or 'yes'  highpass filter
%   cfg.preproc.bpfilter      = 'no' or 'yes'  bandpass filter
%   cfg.preproc.lnfilter      = 'no' or 'yes'  line noise removal using notch filter
%   cfg.preproc.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   cfg.preproc.medianfilter  = 'no' or 'yes'  jump preserving median filter
%   cfg.preproc.lpfreq        = lowpass  frequency in Hz
%   cfg.preproc.hpfreq        = highpass frequency in Hz
%   cfg.preproc.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.preproc.lnfreq        = line noise frequency in Hz, default 50Hz
%   cfg.preproc.lpfiltord     = lowpass  filter order
%   cfg.preproc.hpfiltord     = highpass filter order
%   cfg.preproc.bpfiltord     = bandpass filter order
%   cfg.preproc.lnfiltord     = line noise notch filter order
%   cfg.preproc.medianfiltord = length of median filter
%   cfg.preproc.lpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.preproc.hpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.preproc.bpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.preproc.lpfiltdir     = filter direction, 'twopass' (default) or 'onepass'
%   cfg.preproc.hpfiltdir     = filter direction, 'twopass' (default) or 'onepass'
%   cfg.preproc.bpfiltdir     = filter direction, 'twopass' (default) or 'onepass'
%   cfg.preproc.detrend       = 'no' or 'yes'
%   cfg.preproc.blc           = 'no' or 'yes'
%   cfg.preproc.blcwindow     = [begin end] in seconds, the default is the complete trial
%   cfg.preproc.hilbert       = 'no' or 'yes'
%   cfg.preproc.rectify       = 'no' or 'yes'

% Copyright (C) 2005-2009, Robert Oostenveld
%
% $Log: spikedownsample.m,v $
% Revision 1.46  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.45  2009/01/16 17:21:20  sashae
% added config tracking
%
% Revision 1.44  2009/01/07 10:45:35  roboos
% improved the documentation
%
% Revision 1.43  2009/01/07 10:37:19  roboos
% replaced ADtoUV with cfg.calibration, added explaination to the help
% added some configuration checking
%
% Revision 1.42  2008/12/15 15:07:42  roboos
% give error if cfg.ADtoUV is missing in case of neuralynx_dma or neuralynx_sdma
%
% Revision 1.41  2008/12/15 14:57:06  roboos
% added option cfg.ADtoUV, needed for neuralynx_dma and neuralynx_sdma input datasets
%
% Revision 1.40  2008/10/02 15:32:21  sashae
% replaced call to createsubcfg with checkconfig
%
% Revision 1.39  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.38  2008/09/22 19:43:36  roboos
% switched from write_fcdc_data to write_data
%
% Revision 1.37  2008/07/21 20:10:47  roboos
% updated documentation
%
% Revision 1.36  2008/05/13 15:37:24  roboos
% switched to using read_data/header instead of the read_fcdc_data/header wrapper functions
%
% Revision 1.35  2008/01/30 10:50:41  roboos
% added cfg.timestampdefinition
%
% Revision 1.34  2008/01/14 21:38:01  roboos
% the functionality of the downsampling should have remained the same, but the cfg is not backward compatible
% removed all code related to spikes and mua, only downsampling once remains
% moved all cfg.lfp options into cfg directly
% removed automatic output directory generation
% added cfg.channelprefix to distinguish the channels and the output files
% the function will not return a data structure any more, only the cfg
%
% Revision 1.33  2007/03/21 17:27:30  roboos
% updated documentation

fieldtripdefs
cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the general defaults
if ~isfield(cfg, 'dataset'),            cfg.dataset = [];                 end
if ~isfield(cfg, 'output'),             cfg.output = [];                  end
if ~isfield(cfg, 'channel'),            cfg.channel = 'all';              end
if ~isfield(cfg, 'channelprefix'),      cfg.channelprefix = [];           end
if ~isfield(cfg, 'latency'),            cfg.latency = [0 inf];            end
% set the specific defaults
if ~isfield(cfg, 'fsample'),            cfg.fsample = 1000;               end
%if ~isfield(cfg, 'method'),            cfg.method = [];                  end

if ~isfield(cfg, 'calibration')
  if filetype(cfg.dataset, 'neuralynx_dma') || filetype(cfg.dataset, 'neuralynx_sdma')
    error('You must specify cfg.calibration in case of neuralynx_dma or neuralynx_sdma');
  else
    cfg.calibration = 1;
  end
end

% check that the input cfg is valid for this function
cfg = checkconfig(cfg, 'required', {'method', 'calibration', 'dataformat'});
cfg = checkconfig(cfg, 'forbidden', {'ADtoUV'});

% ensure that the preproc specific options are located in the preproc substructure
cfg = checkconfig(cfg, 'createsubcfg',  {'preproc'});

status = mkdir(cfg.output);
if ~status
  error(sprintf('error creating LFP output dataset %s', cfg.output));
end

% read the header of the completete dataset
hdr = read_header(cfg.dataset);
cfg.channel = channelselection(cfg.channel, hdr.label);
chansel = match_str(hdr.label, cfg.channel);

if strcmp(cfg.timestampdefinition, 'sample')
  % the default would be to keep the original definition of timestamps as determined from looking at the file
  % here the definition of timestamps is changed to correspond with samples at the original sampling rate
  hdr.TimeStampPerSample = 1;
  hdr.FirstTimeStamp     = 1;
  hdr.LastTimeStamp      = hdr.nSamples*hdr.nTrials;
end

if hdr.nSamples<1
  error('the input dataset contains no samples');
elseif length(chansel)<1
  error('the input selection contains no channels');
end

% give some feedback, based on the complete data
fprintf('data contains %10d channels\n', hdr.nChans);
fprintf('selected      %10d channels\n', length(chansel));
numsample = [];
numsegment = size(cfg.latency,1);
for j=1:numsegment
  begsample(j) = max(round(cfg.latency(j,1) * hdr.Fs + 1), 1);
  endsample(j) = min(round(cfg.latency(j,2) * hdr.Fs    ), hdr.nSamples);
  numsample(j) = endsample(j) - begsample(j) + 1;
  cfg.latency(j,1) = (begsample(j)-1)/hdr.Fs;
  cfg.latency(j,2) = (endsample(j)  )/hdr.Fs;
end
numsample = sum(numsample);
fprintf('data contains %10d samples\n', hdr.nSamples);
fprintf('selected      %10d samples in %d segments\n', numsample, numsegment);

s = floor(hdr.nSamples ./ hdr.Fs);
m = floor(s/60);
h = floor(m/60);
m = m - 60*h;
s = s - 60*m - 60*60*h;
fprintf('duration of data      %02dh:%02dm:%02ds\n', h, m, s);

s = floor(numsample ./ hdr.Fs);
m = floor(s/60);
h = floor(m/60);
m = m - 60*h;
s = s - 60*m - 60*60*h;
fprintf('duration of selection %02dh:%02dm:%02ds\n', h, m, s);

fprintf('estimated memory usage %d MB\n', round((numsample*(8+8+2))/(1024^2)));

% determine the exact sampling frequency for the LFP and MUA
fac = hdr.Fs/cfg.fsample;
if strcmp(cfg.method, 'resample');
  % this works with arbitrary changes in the frequency
else
  % this requires an integer downsampling fraction
  fac = round(fac);
  cfg.fsample = hdr.Fs/fac;
end

if numsegment==1
  % determine the timestamps of the first and last sample after downsampling
  tsoff = cast(0, class(hdr.FirstTimeStamp));
  tsbeg = double(hdr.FirstTimeStamp - tsoff);
  tsend = double(hdr.LastTimeStamp  - tsoff);
  tsdif = (tsend-tsbeg)/(hdr.nSamples-1);     % should be the same as hdr.TimeStampPerSample
  % this applies to all samples in the datafile, i.e. not on the selected latency window
  resampled_tsbeg = tsbeg + ((fac-1)/2)*tsdif;
  resampled_tsend = tsend - ((fac-1)/2)*tsdif;
  resampled_tsdif = fac * tsdif;
  % this applies to the selected latency window
  selected_tsbeg = double(hdr.FirstTimeStamp - tsoff) + (begsample-1)*tsdif;
  selected_tsend = double(hdr.FirstTimeStamp - tsoff) + (endsample-1)*tsdif;
  resampled_selected_tsbeg = selected_tsbeg + ((fac-1)/2)*tsdif;
  resampled_selected_tsend = selected_tsend - ((fac-1)/2)*tsdif;
  % these should be integers, but will be casted to uint32 or uint64 later
  resampled_selected_tsbeg = round(resampled_selected_tsbeg);
  resampled_selected_tsend = round(resampled_selected_tsend);
  resampled_tsdif ; % this remains a double
else
  warning('multiple data segments were selected, timestamps are invalid');
end

% process each channel separetely
for i=chansel(:)'
  fprintf('reading channel %d, ''%s''\n', i, hdr.label{i});

  % read the data of a single channel and concatenate into one vector
  org = zeros(1,sum(numsample));
  for j=1:numsegment
    buf = read_data(cfg.dataset, 'header', hdr, 'begsample', begsample(j), 'endsample', endsample(j), 'chanindx', i);

    % apply the optional calibration to the data to ensure that the numbers represent uV
    if cfg.calibration~=1
      buf = cfg.calibration * buf;
    end

    if j==1
      begsegment = 1;
      endsegment = numsample(j);
    else
      begsegment = numsample(j-1) + 1;
      endsegment = sum(numsample(1:j));
    end
    % concatenate the data into one large vector
    org(begsegment:endsegment) = buf;
    clear buf;
  end
  label = hdr.label(i); % this should be cell-array for preproc

  % apply preprocessing and downsample
  fprintf('preprocessing\n');
  dat = preproc(org, label, hdr.Fs, cfg.preproc);
  dat = myresample(dat, hdr.Fs, cfg.fsample, cfg.method);

  chanhdr             = [];
  chanhdr.nChans      = 1;
  chanhdr.nTrials     = 1;             % since continuous
  chanhdr.nSamplesPre = 1;             % since continuous
  chanhdr.Fs          = cfg.fsample;
  chanhdr.nSamples    = length(dat);
  if isempty(cfg.channelprefix)
    % the label should be a cell-array of length one
    chanhdr.label     = hdr.label(i);
  else
    % add a prefix to the channel name
    chanhdr.label     = {[cfg.channelprefix '_' hdr.label{i}]};
  end
  if numsegment==1
    chanhdr.FirstTimeStamp     = resampled_selected_tsbeg;
    chanhdr.LastTimeStamp      = resampled_selected_tsend;
    chanhdr.TimeStampPerSample = resampled_tsdif;
  else
    % these are not defined in case of multiple segments
    warning('multiple data segments were selected, timestamps are invalid');
    chanhdr.FirstTimeStamp     = nan;
    chanhdr.LastTimeStamp      = nan;
    chanhdr.TimeStampPerSample = nan;
  end

  % the output file contains the new channel name
  datafile = fullfile(cfg.output, chanhdr.label{1});  % this is without filename extension
  fprintf('writing to file ''%s''\n', datafile);
  write_data(datafile, dat, 'header', chanhdr, 'dataformat', cfg.dataformat);

  % keep memory tidy
  clear dat

end % for each file

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: spikedownsample.m,v 1.46 2009/01/20 13:01:31 sashae Exp $';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for the resampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dat = myresample(dat, Fold, Fnew, method);
if Fold==Fnew
  return
end
switch method
  case 'resample'
    if ~isa(dat, 'double')
      typ = class(dat);
      dat = typecast(dat, 'double');
      dat = resample(dat, Fnew, Fold);     % this requires a double array
      dat = typecast(dat, typ);
    else
      dat = resample(dat, Fnew, Fold);
    end
  case 'decimate'
    fac = round(Fold/Fnew);
    if ~isa(dat, 'double')
      typ = class(dat);
      dat = typecast(dat, 'double');
      dat = decimate(dat, fac);    % this requires a double array
      dat = typecast(dat, typ);
    else
      dat = decimate(dat, fac);    % this requires a double array
    end
  case 'downsample'
    fac = Fold/Fnew;
    dat = decimate(dat, fac);
  otherwise
    error('unsupported resampling method');
end

