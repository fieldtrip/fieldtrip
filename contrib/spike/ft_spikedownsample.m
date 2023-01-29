function [cfg] = ft_spikedownsample(cfg)

% FT_SPIKEDOWNSAMPLE takes electrophysiological data that was continuoudly
% sampled at 32KHz and preprocesses and downsamples it to obtain the LFP
% data, which can subsequently be processed in more detail.
%
% Use as
%   [cfg] = ft_spikedownsample(cfg)
%
% The configuration should contain
%   cfg.dataset             = string with the input dataset
%   cfg.output              = string with the output dataset (default is determined automatic)
%   cfg.dataformat          = string with the output dataset format, see WRITE_DATA
%   cfg.channel             = Nx1 cell-array with selection of channels (default = 'all'),
%                             see FT_CHANNELSELECTION for details
%   cfg.fsample             = desired sampling frequency in Hz (default = 1000)
%   cfg.method              = resampling method, can be 'resample', 'decimate' or 'downsample'
%   cfg.timestampdefinition = 'orig' or 'sample'
%   cfg.channelprefix       = string, will be added to channel name, e.g. 'lfp' -> 'lfp_ncs001' (default = [])
%   cfg.calibration         = optional scaling factor to apply to the data to convert it in uV, see below
%
% The Neuralynx acquisition system at the FCDC in Nijmegen makes use of
% Plexon headstages which have a amplification of 20x. The data that is
% written by the Neuralynx acquisition software therefore is 20x larger
% than the true microvolt values. When operating FT_SPIKEDOWNSAMPLE on the
% *.ncs files that are recorded with the Neuralynx Cheetah software, the
% calibration should be set to 1/20. The raw dma file (*.nrd) and the
% splitted DMA files contains AD values that are not scaled in uV and
% require an additional factor of 64x. If you operate FT_SPIKEDOWNSAMPLE  on
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
%   cfg.preproc.demean        = 'no' or 'yes'
%   cfg.preproc.baselinewindow = [begin end] in seconds, the default is the complete trial
%   cfg.preproc.hilbert       = 'no' or 'yes'
%   cfg.preproc.rectify       = 'no' or 'yes'

% Copyright (C) 2005-2010, Robert Oostenveld
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
ft_preamble provenance


% set the general defaults
if ~isfield(cfg, 'dataset'),            cfg.dataset = [];                 end
if ~isfield(cfg, 'output'),             cfg.output = [];                  end
if ~isfield(cfg, 'channel'),            cfg.channel = 'all';              end
if ~isfield(cfg, 'channelprefix'),      cfg.channelprefix = [];           end
if ~isfield(cfg, 'latency'),            cfg.latency = [0 inf];            end
if ~isfield(cfg, 'headerformat'),       cfg.headerformat = [];            end 
if ~isfield(cfg, 'dataformat'),         cfg.dataformat = [];              end
% set the specific defaults
if ~isfield(cfg, 'fsample'),            cfg.fsample = 1000;               end
%if ~isfield(cfg, 'method'),            cfg.method = [];                  end
if ~isfield(cfg, 'precision'),          cfg.precision = 'double';         end

if ~isfield(cfg, 'calibration')
  if ft_filetype(cfg.dataset, 'neuralynx_dma') || ft_filetype(cfg.dataset, 'neuralynx_sdma')
    error('You must specify cfg.calibration in case of neuralynx_dma or neuralynx_sdma');
  else
    cfg.calibration = 1;
  end
end

% check that the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', {'method', 'calibration', 'dataformat'});
cfg = ft_checkconfig(cfg, 'forbidden', {'ADtoUV'});

% ensure that the preproc specific options are located in the preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});
cfg.preproc = ft_checkconfig(cfg.preproc, 'renamed', {'blc', 'demean'});
cfg.preproc = ft_checkconfig(cfg.preproc, 'renamed', {'blcwindow', 'baselinewindow'});

status = mkdir(cfg.output);
if ~status
  error('error creating LFP output dataset %s', cfg.output);
end

% read the header of the completete dataset
hdr = ft_read_header(cfg.dataset, 'headerformat', cfg.headerformat);
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
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
    buf = ft_read_data(cfg.dataset, 'header', hdr, 'begsample', begsample(j), 'endsample', endsample(j), 'chanindx', i);

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
  dat = preproc(org, label, offset2time(0, hdr.Fs, size(org,2)), cfg.preproc);
  dat = ft_preproc_resample(dat, hdr.Fs, cfg.fsample, cfg.method);

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
  if strcmp(cfg.dataformat, 'fcdc_matbin')
    chanhdr.precision = cfg.precision;
  end

  % the output file contains the new channel name
  datafile = fullfile(cfg.output, chanhdr.label{1});  % this is without filename extension
  fprintf('writing to file ''%s''\n', datafile);
  ft_write_data(datafile, dat, 'header', chanhdr, 'dataformat', cfg.dataformat);

  % keep memory tidy
  clear dat

end % for each file

% do the general cleanup and bookkeeping at the end of the function

ft_postamble provenance

