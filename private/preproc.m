function [dat, label, time, cfg] = preproc(dat, label, time, cfg, begpadding, endpadding)

% PREPROC applies various preprocessing steps on a piece of EEG/MEG data
% that already has been read from a data file.
%
% This function can serve as a subfunction for all FieldTrip modules that
% want to preprocess the data, such as PREPROCESSING, ARTIFACT_XXX,
% TIMELOCKANALYSIS, etc. It ensures consistent handling of both MEG and EEG
% data and consistency in the use of all preprocessing configuration
% options.
%
% Use as
%   [dat, label, time, cfg] = preproc(dat, label, time, cfg, begpadding, endpadding)
%
% The required input arguments are
%   dat         Nchan x Ntime data matrix
%   label       Nchan x 1 cell-array with channel labels
%   time        Ntime x 1 vector with the latency in seconds
%   cfg         configuration structure, see below
% and the optional input arguments are
%   begpadding  number of samples that was used for padding (see below)
%   endpadding  number of samples that was used for padding (see below)
%
% The output is
%   dat         Nchan x Ntime data matrix
%   label       Nchan x 1 cell-array with channel labels
%   time        Ntime x 1 vector with the latency in seconds
%   cfg         configuration structure, optionally with extra defaults set
%
% Note that the number of input channels and the number of output channels
% can be different, for example when the user specifies that he/she wants
% to add the implicit EEG reference channel to the data matrix.
%
% The filtering of the data can introduce artifacts at the edges, hence it
% is better to pad the data with some extra signal at the begin and end.
% After filtering, this padding is removed and the other preprocessing
% steps are applied to the remainder of the data. The input fields
% begpadding and endpadding should be specified in samples. You can also
% leave them empty, which implies that the data is not padded.
%
% The configuration can contain
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
%   cfg.lpfiltord     = lowpass  filter order (default set in low-level function)
%   cfg.hpfiltord     = highpass filter order (default set in low-level function)
%   cfg.bpfiltord     = bandpass filter order (default set in low-level function)
%   cfg.bsfiltord     = bandstop filter order (default set in low-level function)
%   cfg.medianfiltord = length of median filter
%   cfg.lpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.hpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.bpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.bsfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
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
%   cfg.demean        = 'no' or 'yes'
%   cfg.baselinewindow = [begin end] in seconds, the default is the complete trial
%   cfg.detrend       = 'no' or 'yes', this is done on the complete trial
%   cfg.polyremoval   = 'no' or 'yes', this is done on the complete trial
%   cfg.polyorder     = polynome order (default = 2)
%   cfg.derivative    = 'no' (default) or 'yes', computes the first order derivative of the data
%   cfg.hilbert       = 'no', 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag' or 'angle' (default = 'no')
%   cfg.rectify       = 'no' or 'yes'
%   cfg.precision     = 'single' or 'double' (default = 'double')
%   cfg.absdiff       = 'no' or 'yes', computes absolute derivative (i.e.first derivative then rectify)
%
% Preprocessing options that you should only use for EEG data are
%   cfg.reref         = 'no' or 'yes' (default = 'no')
%   cfg.refchannel    = cell-array with new EEG reference channel(s)
%   cfg.refmethod     = 'avg', 'median', or 'bipolar' (default = 'avg')
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%   cfg.montage       = 'no' or a montage structure (default = 'no')
%
% See also FT_READ_DATA, FT_READ_HEADER

% TODO implement decimation and/or resampling

% Copyright (C) 2004-2012, Robert Oostenveld
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

% compute fsample
fsample = 1./nanmean(diff(time));

if nargin<5 || isempty(begpadding)
  begpadding = 0;
end
if nargin<6 || isempty(endpadding)
  endpadding = 0;
end

if iscell(cfg)
  % recurse over the subsequent preprocessing stages
  if begpadding>0 || endpadding>0
    ft_error('multiple preprocessing stages are not supported in combination with filter padding');
  end
  for i=1:length(cfg)
    tmpcfg = cfg{i};
    if nargout==1
      [dat                     ] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
    elseif nargout==2
      [dat, label              ] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
    elseif nargout==3
      [dat, label, time        ] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
    elseif nargout==4
      [dat, label, time, tmpcfg] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
      cfg{i} = tmpcfg;
    end
  end
  % ready with recursing over the subsequent preprocessing stages
  return
end

% set the defaults for the rereferencing options
cfg.reref =                ft_getopt(cfg, 'reref', 'no');
cfg.refchannel =           ft_getopt(cfg, 'refchannel', {});
cfg.refmethod =            ft_getopt(cfg, 'refmethod', 'avg');
cfg.implicitref =          ft_getopt(cfg, 'implicitref', []);
% set the defaults for the signal processing options
cfg.polyremoval =          ft_getopt(cfg, 'polyremoval', 'no');
cfg.polyorder =            ft_getopt(cfg, 'polyorder', 2);
cfg.detrend =              ft_getopt(cfg, 'detrend', 'no');
cfg.demean =               ft_getopt(cfg, 'demean', 'no');
cfg.baselinewindow =       ft_getopt(cfg, 'baselinewindow', 'all');
cfg.dftfilter =            ft_getopt(cfg, 'dftfilter', 'no');
cfg.lpfilter =             ft_getopt(cfg, 'lpfilter', 'no');
cfg.hpfilter =             ft_getopt(cfg, 'hpfilter', 'no');
cfg.bpfilter =             ft_getopt(cfg, 'bpfilter', 'no');
cfg.bsfilter =             ft_getopt(cfg, 'bsfilter', 'no');
cfg.lpfiltord =            ft_getopt(cfg, 'lpfiltord', []);
cfg.hpfiltord =            ft_getopt(cfg, 'hpfiltord', []);
cfg.bpfiltord =            ft_getopt(cfg, 'bpfiltord', []);
cfg.bsfiltord =            ft_getopt(cfg, 'bsfiltord', []);
cfg.lpfilttype =           ft_getopt(cfg, 'lpfilttype', 'but');
cfg.hpfilttype =           ft_getopt(cfg, 'hpfilttype', 'but');
cfg.bpfilttype =           ft_getopt(cfg, 'bpfilttype', 'but');
cfg.bsfilttype =           ft_getopt(cfg, 'bsfilttype', 'but');
if strcmp(cfg.lpfilttype, 'firws'), cfg.lpfiltdir = ft_getopt(cfg, 'lpfiltdir', 'onepass-zerophase'); else, cfg.lpfiltdir = ft_getopt(cfg, 'lpfiltdir', 'twopass'); end
if strcmp(cfg.hpfilttype, 'firws'), cfg.hpfiltdir = ft_getopt(cfg, 'hpfiltdir', 'onepass-zerophase'); else, cfg.hpfiltdir = ft_getopt(cfg, 'hpfiltdir', 'twopass'); end
if strcmp(cfg.bpfilttype, 'firws'), cfg.bpfiltdir = ft_getopt(cfg, 'bpfiltdir', 'onepass-zerophase'); else, cfg.bpfiltdir = ft_getopt(cfg, 'bpfiltdir', 'twopass'); end
if strcmp(cfg.bsfilttype, 'firws'), cfg.bsfiltdir = ft_getopt(cfg, 'bsfiltdir', 'onepass-zerophase'); else, cfg.bsfiltdir = ft_getopt(cfg, 'bsfiltdir', 'twopass'); end
cfg.lpinstabilityfix =     ft_getopt(cfg, 'lpinstabilityfix', 'no');
cfg.hpinstabilityfix =     ft_getopt(cfg, 'hpinstabilityfix', 'no');
cfg.bpinstabilityfix =     ft_getopt(cfg, 'bpinstabilityfix', 'no');
cfg.bsinstabilityfix =     ft_getopt(cfg, 'bsinstabilityfix', 'no');
cfg.lpfiltdf =             ft_getopt(cfg, 'lpfiltdf',  []);
cfg.hpfiltdf =             ft_getopt(cfg, 'hpfiltdf', []);
cfg.bpfiltdf =             ft_getopt(cfg, 'bpfiltdf', []);
cfg.bsfiltdf =             ft_getopt(cfg, 'bsfiltdf', []);
cfg.lpfiltwintype =        ft_getopt(cfg, 'lpfiltwintype', 'hamming');
cfg.hpfiltwintype =        ft_getopt(cfg, 'hpfiltwintype', 'hamming');
cfg.bpfiltwintype =        ft_getopt(cfg, 'bpfiltwintype', 'hamming');
cfg.bsfiltwintype =        ft_getopt(cfg, 'bsfiltwintype', 'hamming');
cfg.lpfiltdev =            ft_getopt(cfg, 'lpfiltdev', []);
cfg.hpfiltdev =            ft_getopt(cfg, 'hpfiltdev', []);
cfg.bpfiltdev =            ft_getopt(cfg, 'bpfiltdev', []);
cfg.bsfiltdev =            ft_getopt(cfg, 'bsfiltdev', []);
cfg.plotfiltresp =         ft_getopt(cfg, 'plotfiltresp', 'no');
cfg.usefftfilt =           ft_getopt(cfg, 'usefftfilt', 'no');
cfg.medianfilter =         ft_getopt(cfg, 'medianfilter ', 'no');
cfg.medianfiltord =        ft_getopt(cfg, 'medianfiltord', 9);
cfg.dftfreq =              ft_getopt(cfg, 'dftfreq', [50 100 150]);
cfg.hilbert =              ft_getopt(cfg, 'hilbert', 'no');
cfg.derivative =           ft_getopt(cfg, 'derivative', 'no');
cfg.rectify =              ft_getopt(cfg, 'rectify', 'no');
cfg.boxcar =               ft_getopt(cfg, 'boxcar', 'no');
cfg.absdiff =              ft_getopt(cfg, 'absdiff', 'no');
cfg.precision =            ft_getopt(cfg, 'precision', []);
cfg.conv =                 ft_getopt(cfg, 'conv', 'no');
cfg.montage =              ft_getopt(cfg, 'montage', 'no');
cfg.dftinvert =            ft_getopt(cfg, 'dftinvert', 'no');
cfg.standardize =          ft_getopt(cfg, 'standardize', 'no');
cfg.denoise =              ft_getopt(cfg, 'denoise', '');
cfg.subspace =             ft_getopt(cfg, 'subspace', []);
cfg.custom =               ft_getopt(cfg, 'custom', '');
cfg.resample =             ft_getopt(cfg, 'resample', '');

% test whether the MATLAB signal processing toolbox is available
if strcmp(cfg.medianfilter, 'yes') && ~ft_hastoolbox('signal')
  ft_error('median filtering requires the MATLAB signal processing toolbox');
end

% do a sanity check on the filter configuration
if strcmp(cfg.bpfilter, 'yes') && ...
    (strcmp(cfg.hpfilter, 'yes') || strcmp(cfg.lpfilter,'yes'))
  ft_error('you should not apply both a bandpass AND a lowpass/highpass filter');
end

% do a sanity check on the hilbert transform configuration
if strcmp(cfg.hilbert, 'yes') && ~strcmp(cfg.bpfilter, 'yes')
  ft_warning('Hilbert transform should be applied in conjunction with bandpass filter')
end

% do a sanity check on hilbert and rectification
if strcmp(cfg.hilbert, 'yes') && strcmp(cfg.rectify, 'yes')
  ft_error('Hilbert transform and rectification should not be applied both')
end

% do a sanity check on the rereferencing/montage
if ~strcmp(cfg.reref, 'no') && ~strcmp(cfg.montage, 'no')
  ft_error('cfg.reref and cfg.montage are mutually exclusive')
end

% lnfilter is no longer used
if isfield(cfg, 'lnfilter') && strcmp(cfg.lnfilter, 'yes')
  ft_error('line noise filtering using the option cfg.lnfilter is not supported any more, use cfg.bsfilter instead')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the rereferencing in case of EEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.implicitref) && ~any(match_str(cfg.implicitref,label))
  label = {label{:} cfg.implicitref}';
  dat(end+1,:) = 0;
end

if strcmp(cfg.reref, 'yes')
  if strcmp(cfg.refmethod, 'bipolar')
    % this is implemented as a montage that the user does not get to see
    tmpcfg = keepfields(cfg, {'refmethod', 'implicitref', 'refchannel', 'channel'});
    tmpcfg.showcallinfo = 'no';
    montage = ft_prepare_montage(tmpcfg);
    % convert the data temporarily to a raw structure
    tmpdata.trial = {dat};
    tmpdata.time  = {time};
    tmpdata.label = label;
    % apply the montage to the data
    tmpdata = ft_apply_montage(tmpdata, montage, 'feedback', 'none');
    dat   = tmpdata.trial{1}; % the number of channels can have changed
    label = tmpdata.label;    % the output channels can be different than the input channels
    clear tmpdata
  else
    % mean or median based derivation of specified or all channels
    cfg.refchannel = ft_channelselection(cfg.refchannel, label);
    refindx = match_str(label, cfg.refchannel);
    if isempty(refindx)
      ft_error('reference channel was not found')
    end
    dat = ft_preproc_rereference(dat, refindx, cfg.refmethod);
  end
end

if ~strcmp(cfg.montage, 'no') && ~isempty(cfg.montage)
  % convert the data temporarily to a raw structure
  tmpdata.trial = {dat};
  tmpdata.time  = {time};
  tmpdata.label = label;
  % apply the montage to the data
  tmpdata = ft_apply_montage(tmpdata, cfg.montage, 'feedback', 'none');
  dat   = tmpdata.trial{1}; % the number of channels can have changed
  label = tmpdata.label;    % the output channels can be different than the input channels
  clear tmpdata
end

if any(any(isnan(dat)))
  % filtering is not possible for at least a selection of the data
  ft_warning('data contains NaNs, no filtering or preprocessing applied');
  
else
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do the filtering on the padded data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isempty(cfg.denoise)
    hflag    = isfield(cfg.denoise, 'hilbert') && strcmp(cfg.denoise.hilbert, 'yes');
    datlabel = match_str(label, cfg.denoise.channel);
    reflabel = match_str(label, cfg.denoise.refchannel);
    tmpdat   = ft_preproc_denoise(dat(datlabel,:), dat(reflabel,:), hflag);
    dat(datlabel,:) = tmpdat;
  end
  
  % The filtering should in principle be done prior to the demeaning to
  % ensure that the resulting mean over the baseline window will be
  % guaranteed to be zero (even if there are filter artifacts).
  % However, the filtering benefits from the data being pulled towards zero,
  % causing less edge artifacts. That is why we start by removing the slow
  % drift, then filter, and then repeat the demean/detrend/polyremove.
  if strcmp(cfg.polyremoval, 'yes')
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_polyremoval(dat, cfg.polyorder, begsample, endsample); % this will also demean and detrend
  elseif strcmp(cfg.detrend, 'yes')
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_polyremoval(dat, 1, begsample, endsample); % this will also demean
  elseif strcmp(cfg.demean, 'yes')
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_polyremoval(dat, 0, begsample, endsample);
  end
  
  if strcmp(cfg.medianfilter, 'yes'), dat = ft_preproc_medianfilter(dat, cfg.medianfiltord); end
  if strcmp(cfg.lpfilter, 'yes'),     dat = ft_preproc_lowpassfilter(dat, fsample, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir, cfg.lpinstabilityfix, cfg.lpfiltdf, cfg.lpfiltwintype, cfg.lpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
  if strcmp(cfg.hpfilter, 'yes'),     dat = ft_preproc_highpassfilter(dat, fsample, cfg.hpfreq, cfg.hpfiltord, cfg.hpfilttype, cfg.hpfiltdir, cfg.hpinstabilityfix, cfg.hpfiltdf, cfg.hpfiltwintype, cfg.hpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
  if strcmp(cfg.bpfilter, 'yes'),     dat = ft_preproc_bandpassfilter(dat, fsample, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir, cfg.bpinstabilityfix, cfg.bpfiltdf, cfg.bpfiltwintype, cfg.bpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
  if strcmp(cfg.bsfilter, 'yes')
    for i=1:size(cfg.bsfreq,1)
      % apply a bandstop filter for each of the specified bands, i.e. cfg.bsfreq should be Nx2
      dat = ft_preproc_bandstopfilter(dat, fsample, cfg.bsfreq(i,:), cfg.bsfiltord, cfg.bsfilttype, cfg.bsfiltdir, cfg.bsinstabilityfix, cfg.bsfiltdf, cfg.bsfiltwintype, cfg.bsfiltdev, cfg.plotfiltresp, cfg.usefftfilt);
    end
  end
  if strcmp(cfg.polyremoval, 'yes')
    % the begin and endsample of the polyremoval period correspond to the complete data minus padding
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_polyremoval(dat, cfg.polyorder, begsample, endsample);
  end
  if strcmp(cfg.detrend, 'yes')
    % the begin and endsample of the detrend period correspond to the complete data minus padding
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_detrend(dat, begsample, endsample);
  end
  if strcmp(cfg.demean, 'yes')
    if ischar(cfg.baselinewindow) && strcmp(cfg.baselinewindow, 'all')
      % the begin and endsample of the baseline period correspond to the complete data minus padding
      nsamples  = size(dat,2);
      begsample = 1        + begpadding;
      endsample = nsamples - endpadding;
      dat       = ft_preproc_baselinecorrect(dat, begsample, endsample);
    else
      % determine the begin and endsample of the baseline period and baseline correct for it
      begsample = nearest(time, cfg.baselinewindow(1));
      endsample = nearest(time, cfg.baselinewindow(2));
      dat       = ft_preproc_baselinecorrect(dat, begsample, endsample);
    end
  end
  if strcmp(cfg.dftfilter, 'yes')
    datorig = dat;
    optarg = {};
    if isfield(cfg, 'dftreplace')
      optarg = cat(2, optarg, {'dftreplace', cfg.dftreplace});
      if strcmp(cfg.dftreplace, 'neighbour') && (begpadding>0 || endpadding>0)
        ft_error('Padding by data mirroring is not supported for spectrum interpolation.');
      end
    end
    if isfield(cfg, 'dftbandwidth')
      optarg = cat(2, optarg, {'dftbandwidth', cfg.dftbandwidth});
    end
    if isfield(cfg, 'dftneighbourwidth')
      optarg = cat(2, optarg, {'dftneighbourwidth', cfg.dftneighbourwidth});
    end
    dat     = ft_preproc_dftfilter(dat, fsample, cfg.dftfreq, optarg{:});
    if strcmp(cfg.dftinvert, 'yes')
      dat = datorig - dat;
    end
  end
  if ~strcmp(cfg.hilbert, 'no')
    dat = ft_preproc_hilbert(dat, cfg.hilbert);
  end
  if strcmp(cfg.rectify, 'yes')
    dat = ft_preproc_rectify(dat);
  end
  if isnumeric(cfg.boxcar)
    numsmp = round(cfg.boxcar*fsample);
    if ~rem(numsmp,2)
      % the kernel should have an odd number of samples
      numsmp = numsmp+1;
    end
    % kernel = ones(1,numsmp) ./ numsmp;
    % dat    = convn(dat, kernel, 'same');
    dat = ft_preproc_smooth(dat, numsmp); % better edge behavior
  end
  if isnumeric(cfg.conv)
    kernel = (cfg.conv(:)'./sum(cfg.conv));
    if ~rem(length(kernel),2)
      kernel = [kernel 0];
    end
    dat = convn(dat, kernel, 'same');
  end
  if strcmp(cfg.derivative, 'yes')
    dat = ft_preproc_derivative(dat, 1);
  end
  if strcmp(cfg.absdiff, 'yes')
    % this implements abs(diff(data), which is required for jump detection
    dat = abs([diff(dat, 1, 2) zeros(size(dat,1),1)]);
  end
  if strcmp(cfg.standardize, 'yes')
    dat = ft_preproc_standardize(dat, 1, size(dat,2));
  end
  if ~isempty(cfg.subspace)
    dat = ft_preproc_subspace(dat, cfg.subspace);
  end
  if ~isempty(cfg.custom)
    if ~isfield(cfg.custom, 'nargout')
      cfg.custom.nargout = 1;
    end
    if cfg.custom.nargout==1
      dat = feval(cfg.custom.funhandle, dat, cfg.custom.varargin);
    elseif cfg.custom.nargout==2
      [dat, time] = feval(cfg.custom.funhandle, dat, cfg.custom.varargin);
    end
  end
  if strcmp(cfg.resample, 'yes')
    if ~isfield(cfg, 'resamplefs')
      cfg.resamplefs = fsample./2;
    end
    if ~isfield(cfg, 'resamplemethod')
      cfg.resamplemethod = 'resample';
    end
    [dat               ] = ft_preproc_resample(dat,  fsample, cfg.resamplefs, cfg.resamplemethod);
    [time, dum, fsample] = ft_preproc_resample(time, fsample, cfg.resamplefs, cfg.resamplemethod);
  end
  if ~isempty(cfg.precision)
    % convert the data to another numeric precision, i.e. double, single or int32
    dat = cast(dat, cfg.precision);
  end
end % if any(isnan)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the filter padding and do the preprocessing on the remaining trial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if begpadding~=0 || endpadding~=0
  dat = ft_preproc_padding(dat, 'remove', begpadding, endpadding);
  if strcmp(cfg.demean, 'yes') || nargout>2
    time = ft_preproc_padding(time, 'remove', begpadding, endpadding);
  end
end
