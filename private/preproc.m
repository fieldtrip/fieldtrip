function [dat, label, time, cfg] = preproc(dat, label, fsample, cfg, offset, begpadding, endpadding);

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
%   [dat, label, time, cfg] = preproc(dat, label, fsample, cfg, offset, begpadding, endpadding)
%
% The required input arguments are
%   dat         Nchan x Ntime data matrix
%   label       Nchan x 1 cell-array with channel labels
%   cfg         configuration structure, see below
%   fsample     sampling frequency
% and the optional input arguments are
%   offset      of the first datasample (see below)
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
% The offset field specifies the difference in the latency of the beginning
% of the data relative to the occurence of the trigger (expressed in
% samples). An offset of 0 means that the first sample of the trial
% corresponds with the trigger. A positive offset indicates that the first
% sample is later than the triger, a negative offset indicates that the
% trial begins before the trigger. The offset should be specified EXCLUDING
% the filter padding at the begin of the data. You can leave it empty,
% which implies that the first sample of the data corresponds with latency 0.
%
% The filtering of the data can introduce artifacts at the edges, hence it
% is better to pad the data with some extra signal at the begin and end.
% After filtering, this padding is removed and the other preprocessing
% steps are applied to the remainder of the data. The input fields
% begpadding and endpadding should be specified in samples. You can also
% leave them empty, which implies that the data is not padding.
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
%   cfg.lpfiltord     = lowpass  filter order
%   cfg.hpfiltord     = highpass filter order
%   cfg.bpfiltord     = bandpass filter order
%   cfg.bsfiltord     = bandstop filter order
%   cfg.medianfiltord = length of median filter
%   cfg.lpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.hpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.bpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.bsfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.lpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.hpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.bpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.bsfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
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
% See also READ_DATA, READ_HEADER

% TODO implement decimation and/or resampling

% Copyright (C) 2004-2009, Robert Oostenveld
%
% $Log: preproc.m,v $
% Revision 1.38  2009/09/30 12:58:53  jansch
% added optional call to preproc_denoise and preproc_subspace
%
% Revision 1.37  2009/08/05 08:36:09  roboos
% don't fill up the complete data with nans if a nan is detected
% the behaviour of not filtering and giving a warning remains the same
%
% Revision 1.36  2009/06/17 10:12:26  roboos
% cfg.montage=[] should also work (just like 'no')
%
% Revision 1.35  2009/03/13 13:24:00  jansch
% added support for preproc_denoise
%
% Revision 1.34  2009/03/11 11:26:43  roboos
% updated documentation and copyrights
%
% Revision 1.33  2008/10/10 09:54:53  jansch
% added (undocumented) option dftinvert which results in the dftfilter being
% a very sharp bandpass instead of a notch filter. added (undocumented) option
% conv. fixed small bug in derivative
%
% Revision 1.32  2008/07/08 08:05:44  sashae
% fixed small bug in derivative
%
% Revision 1.31  2008/07/07 14:59:14  sashae
% now using preproc_modules for low-level preprocessing functions
% lnfilter is no longer supported
% updated documentation
%
% Revision 1.30  2008/06/26 07:55:38  roboos
% apply_montage needs cell+label structure and not matrix (thanks to Conrado)
%
% Revision 1.29  2008/06/25 06:37:17  roboos
% change in whitespace
%
% Revision 1.28  2008/06/24 12:47:12  roboos
% added cfg.montage as alternative for rereferencing
%
% Revision 1.27  2008/04/09 14:12:41  roboos
% test isnan over both dimensions of data
%
% Revision 1.26  2007/11/27 16:42:17  roboos
% unwrap angle for hilbert, fixed typo for absimag
%
% Revision 1.25  2007/11/14 13:16:14  roboos
% loop over multiple frequency bands for bandstopfilter, thanks to Saskia
%
% Revision 1.24  2007/11/08 10:09:50  roboos
% allow other versions of the hilbert transformed signal to be computed (e.g. complex, real, imag)
%
% Revision 1.23  2007/10/16 12:38:10  roboos
% do not process the data if there are NaNs, give warning and replace the data with all NaNs
%
% Revision 1.22  2007/09/20 12:59:40  roboos
% updated documentation, added cfg.precision for typecasting data to double or single
%
% Revision 1.21  2007/08/01 15:10:47  ingnie
% added option polyremoval, renamed internal variable blcbeg/endsample into
% beg/endsample
%
% Revision 1.20  2007/04/16 19:04:11  roboos
% removed an excess "end", thanks to Jo
%
% Revision 1.19  2007/04/16 16:10:36  roboos
% loop over all frequencies specified in the cfg.lnfreq (for notch filtering, 2Hz wide)
% added support bandstop filtering (cfg.bsfilter) for when the user wants to specify a wider stop-band
%
% Revision 1.18  2007/01/17 17:05:10  roboos
% use hastoolbox('signal')
%
% Revision 1.17  2006/08/31 07:56:05  roboos
% added onepass-reverse filter to documentation
%
% Revision 1.16  2006/06/14 12:45:58  roboos
% removed documentation of non-functional option cfg.lnfilttype
% added support for onepass and twopass filtering using cfg.lpfiltdir etc.
%
% Revision 1.15  2006/06/14 11:56:30  roboos
% added support for multiple preprocessing stages, achieved by using a cell-array as input (each cell containing a seperate cfg)
%
% Revision 1.14  2006/04/25 20:20:50  roboos
% moved some of the sanity checks from preprocessing to private/preproc
% reinserted the default of some of the cfg settings, since that was broken
%
% Revision 1.13  2006/02/06 20:41:35  roboos
% fixed bug: padding should be removed from time axis
% fixed bug: blc window selection should be done on time axis keeping padding in mind
%
% Revision 1.12  2006/01/18 15:00:00  jansch
% moved the removal of the filter-padding to the end. replaced conv and for-loop
% by convn. implemented mydetrend as a subfunction, which behaves similar to blc.
%
% Revision 1.11  2006/01/17 14:07:43  roboos
% moved rereferencing for EEG all the way to the top, prior to filtering
% implemented cfg.absdiff=yes|no, which does abs(diff(data)) to ensure the order of these two operations (important for jump detection)
%
% Revision 1.10  2006/01/12 16:02:27  roboos
% fixed bug in boxcar, loop over conv() for multiple channels
%
% Revision 1.9  2005/12/20 13:21:46  roboos
% added boxcar convolution as proccessing method
% added temporal derivative as processing method
%
% Revision 1.8  2005/12/02 08:58:29  roboos
% made construction of the time axis optional (can take large amount of memory)
%
% Revision 1.7  2005/11/23 10:44:14  roboos
% added bp/lp/hp/lnfilttype to the documentation
%
% Revision 1.6  2005/09/08 16:56:00  roboos
% only remove padding if unequal to zero
% changed location in file where the time axis is computed
%
% Revision 1.5  2005/09/02 13:17:50  roboos
% Changed dftfilter, loop over all specified frequencies instead of explicitely taking the 2x and 3x harmonics
% Changed default for dftfilter, it now used cfg.dftfreq instead of cfg.lnfreq
% The default is cfg.dftfreq=[50 100 150]
%
% Revision 1.4  2005/05/17 17:50:49  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.3  2005/05/02 08:17:46  roboos
% implemented suggestion of Christian Forkstam: only add implicit reference channel if it is not yet present in the data
%
% Revision 1.2  2005/01/21 09:53:11  roboos
% implemented median filter in preproc, updated help
%
% Revision 1.1  2004/12/09 17:22:28  roboos
% initial version, replicates all preprocessing steps from the large preprocessing function such as filtering, detrending and re-referencing.
%

if nargin<5 || isempty(offset)
  offset = 0;
end
if nargin<6 || isempty(begpadding)
  begpadding = 0;
end
if nargin<7 || isempty(endpadding)
  endpadding = 0;
end

if iscell(cfg)
  % recurse over the subsequent preprocessing stages
  if begpadding>0 || endpadding>0
    error('multiple preprocessing stages are not supported in combination with filter padding');
  end
  for i=1:length(cfg)
    tmpcfg = cfg{i};
    if nargout==1
      [dat                     ] = preproc(dat, label, fsample, tmpcfg, offset, begpadding, endpadding);
    elseif nargout==2
      [dat, label              ] = preproc(dat, label, fsample, tmpcfg, offset, begpadding, endpadding);
    elseif nargout==3
      [dat, label, time        ] = preproc(dat, label, fsample, tmpcfg, offset, begpadding, endpadding);
    elseif nargout==4
      [dat, label, time, tmpcfg] = preproc(dat, label, fsample, tmpcfg, offset, begpadding, endpadding);
      cfg{i} = tmpcfg;
    end
  end
  % ready with recursing over the subsequent preprocessing stages
  return
end

% set the defaults for the rereferencing options
if ~isfield(cfg, 'reref'),        cfg.reref = 'no';             end
if ~isfield(cfg, 'refchannel'),   cfg.refchannel = {};          end
if ~isfield(cfg, 'implicitref'),  cfg.implicitref = [];         end
% set the defaults for the signal processing options
if ~isfield(cfg, 'polyremoval'),  cfg.polyremoval = 'no';       end
if ~isfield(cfg, 'polyorder'),    cfg.polyorder = 2;            end
if ~isfield(cfg, 'detrend'),      cfg.detrend = 'no';           end
if ~isfield(cfg, 'blc'),          cfg.blc = 'no';               end
if ~isfield(cfg, 'blcwindow'),    cfg.blcwindow = 'all';        end
if ~isfield(cfg, 'dftfilter'),    cfg.dftfilter = 'no';         end
if ~isfield(cfg, 'lpfilter'),     cfg.lpfilter = 'no';          end
if ~isfield(cfg, 'hpfilter'),     cfg.hpfilter = 'no';          end
if ~isfield(cfg, 'bpfilter'),     cfg.bpfilter = 'no';          end
if ~isfield(cfg, 'bsfilter'),     cfg.bsfilter = 'no';          end
if ~isfield(cfg, 'lpfiltord'),    cfg.lpfiltord = 6;            end
if ~isfield(cfg, 'hpfiltord'),    cfg.hpfiltord = 6;            end
if ~isfield(cfg, 'bpfiltord'),    cfg.bpfiltord = 4;            end
if ~isfield(cfg, 'bsfiltord'),    cfg.bsfiltord = 4;            end
if ~isfield(cfg, 'lpfilttype'),   cfg.lpfilttype = 'but';       end
if ~isfield(cfg, 'hpfilttype'),   cfg.hpfilttype = 'but';       end
if ~isfield(cfg, 'bpfilttype'),   cfg.bpfilttype = 'but';       end
if ~isfield(cfg, 'bsfilttype'),   cfg.bsfilttype = 'but';       end
if ~isfield(cfg, 'lpfiltdir'),    cfg.lpfiltdir = 'twopass';    end
if ~isfield(cfg, 'hpfiltdir'),    cfg.hpfiltdir = 'twopass';    end
if ~isfield(cfg, 'bpfiltdir'),    cfg.bpfiltdir = 'twopass';    end
if ~isfield(cfg, 'bsfiltdir'),    cfg.bsfiltdir = 'twopass';    end
if ~isfield(cfg, 'medianfilter'), cfg.medianfilter  = 'no';     end
if ~isfield(cfg, 'medianfiltord'),cfg.medianfiltord = 9;        end
if ~isfield(cfg, 'dftfreq')       cfg.dftfreq = [50 100 150];   end
if ~isfield(cfg, 'hilbert'),      cfg.hilbert = 'no';           end
if ~isfield(cfg, 'derivative'),   cfg.derivative = 'no';        end
if ~isfield(cfg, 'rectify'),      cfg.rectify = 'no';           end
if ~isfield(cfg, 'boxcar'),       cfg.boxcar = 'no';            end
if ~isfield(cfg, 'absdiff'),      cfg.absdiff = 'no';           end
if ~isfield(cfg, 'precision'),    cfg.precision = [];           end
if ~isfield(cfg, 'conv'),         cfg.conv = 'no';              end
if ~isfield(cfg, 'montage'),      cfg.montage = 'no';           end
if ~isfield(cfg, 'dftinvert'),    cfg.dftinvert = 'no';         end
if ~isfield(cfg, 'standardize'),  cfg.standardize = 'no';       end
if ~isfield(cfg, 'denoise'),      cfg.denoise = '';             end
if ~isfield(cfg, 'subspace'),     cfg.subspace = [];            end

% test whether the Matlab signal processing toolbox is available
if strcmp(cfg.medianfilter, 'yes') && ~hastoolbox('signal')
  error('median filtering requires the Matlab signal processing toolbox');
end

% do a sanity check on the filter configuration
if strcmp(cfg.bpfilter, 'yes') && ...
    (strcmp(cfg.hpfilter, 'yes') || strcmp(cfg.lpfilter,'yes')),
  error('you should not apply both a bandpass AND a lowpass/highpass filter');
end

% do a sanity check on the hilbert transform configuration
if strcmp(cfg.hilbert, 'yes') && ~strcmp(cfg.bpfilter, 'yes')
  error('hilbert transform should be applied in conjunction with bandpass filter')
end

% do a sanity check on hilbert and rectification
if strcmp(cfg.hilbert, 'yes') && strcmp(cfg.rectify, 'yes')
  error('hilbert transform and rectification should not be applied both')
end

% do a sanity check on the rereferencing/montage
if ~strcmp(cfg.reref, 'no') && ~strcmp(cfg.montage, 'no')
  error('cfg.reref and cfg.montage are mutually exclusive')
end

% lnfilter is no longer used
if isfield(cfg, 'lnfilter') && strcmp(cfg.lnfilter, 'yes')
  error('line noise filtering using the option cfg.lnfilter is not supported any more, use cfg.bsfilter instead')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the rereferencing in case of EEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.implicitref) && ~any(strmatch(cfg.implicitref,label))
  label = {label{:} cfg.implicitref};
  dat(end+1,:) = 0;
end

if strcmp(cfg.reref, 'yes'),
  cfg.refchannel = channelselection(cfg.refchannel, label);
  refindx = match_str(label, cfg.refchannel);
  if isempty(refindx)
    error('reference channel was not found')
  end
  dat = preproc_rereference(dat, refindx);
end

if ~strcmp(cfg.montage, 'no') && ~isempty(cfg.montage)
  % this is an alternative approach for rereferencing, with arbitrary complex linear combinations of channels
  tmp.trial = {dat};
  tmp.label = label;
  tmp = apply_montage(tmp, cfg.montage);
  dat = tmp.trial{1};
  label = tmp.label;
  clear tmp
end

if any(any(isnan(dat)))
  % filtering is not possible for at least a selection of the data
  if nargout>2
    nsamples = size(dat,2);
    time = (offset - begpadding + (0:(nsamples-1)))/fsample;
  end
  warning('data contains NaNs, no filtering applied');
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the filtering on the padded data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.denoise),
  hflag    = isfield(cfg.denoise, 'hilbert') && strcmp(cfg.denoise.hilbert, 'yes');
  datlabel = match_str(label, cfg.denoise.channel);
  reflabel = match_str(label, cfg.denoise.refchannel);
  tmpdat   = preproc_denoise(dat(datlabel,:), dat(reflabel,:), hflag);
  dat(datlabel,:) = tmpdat;
end
if strcmp(cfg.medianfilter, 'yes'), dat = preproc_medianfilter(dat, cfg.medianfiltord); end
if strcmp(cfg.lpfilter, 'yes'),     dat = preproc_lowpassfilter(dat, fsample, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir); end
if strcmp(cfg.hpfilter, 'yes'),     dat = preproc_highpassfilter(dat, fsample, cfg.hpfreq, cfg.hpfiltord, cfg.hpfilttype, cfg.hpfiltdir); end
if strcmp(cfg.bpfilter, 'yes'),     dat = preproc_bandpassfilter(dat, fsample, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir); end
if strcmp(cfg.bsfilter, 'yes')
  for i=1:size(cfg.bsfreq,1)
    % apply a bandstop filter for each of the specified bands, i.e. cfg.bsfreq should be Nx2
    dat = preproc_bandstopfilter(dat, fsample, cfg.bsfreq(i,:), cfg.bsfiltord, cfg.bsfilttype, cfg.bsfiltdir);
  end
end
if strcmp(cfg.dftfilter, 'yes')
  datorig = dat;
  for i=1:length(cfg.dftfreq)
    % filter out the 50Hz noise, optionally also the 100 and 150 Hz harmonics
    dat = preproc_dftfilter(dat, fsample, cfg.dftfreq(i));
  end
  if strcmp(cfg.dftinvert, 'yes'),
    dat = datorig - dat;
  end
end
if strcmp(cfg.polyremoval, 'yes')
  nsamples     = size(dat,2);
  % the begin and endsample of the polyremoval period correspond to the complete data minus padding
  begsample = 1        + begpadding;
  endsample = nsamples - endpadding;
  dat = polyremoval(dat, cfg.polyorder, begsample, endsample);
end
if strcmp(cfg.detrend, 'yes')
  nsamples     = size(dat,2);
  % the begin and endsample of the detrend period correspond to the complete data minus padding
  begsample = 1        + begpadding;
  endsample = nsamples - endpadding;
  dat = preproc_detrend(dat, begsample, endsample);
end
if strcmp(cfg.blc, 'yes') || nargout>2
  % determine the complete time axis for the baseline correction
  % but only construct it when really needed, since it takes up a large amount of memory
  % the time axis should include the filter padding
  nsamples = size(dat,2);
  time = (offset - begpadding + (0:(nsamples-1)))/fsample;
end
if strcmp(cfg.blc, 'yes')
  if isstr(cfg.blcwindow) && strcmp(cfg.blcwindow, 'all')
    % the begin and endsample of the baseline period correspond to the complete data minus padding
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat       = preproc_baselinecorrect(dat, begsample, endsample);
  else
    % determine the begin and endsample of the baseline period and baseline correct for it
    begsample = nearest(time, cfg.blcwindow(1));
    endsample = nearest(time, cfg.blcwindow(2));
    dat       = preproc_baselinecorrect(dat, begsample, endsample);
  end
end
if ~strcmp(cfg.hilbert, 'no')
  dat = preproc_hilbert(dat, cfg.hilbert);
end
if strcmp(cfg.rectify, 'yes'),
  dat = preproc_rectify(dat);
end
if isnumeric(cfg.boxcar)
  numsmp = round(cfg.boxcar*fsample);
  if ~rem(numsmp,2)
    % the kernel should have an odd number of samples
    numsmp = numsmp+1;
  end
  kernel = ones(1,numsmp) ./ numsmp;
  %begsmp = (numsmp-1)/2 + 1;
  %endsmp = (numsmp-1)/2 + size(dat,2);
  %for i=1:size(dat,1)
  %  tmp = conv(dat(i,:), kernel);
  %  % remove the kernel padding at the edges
  %  dat(i,:) = tmp(begsmp:endsmp);
  %end
  dat = convn(dat, kernel, 'same');
end
if isnumeric(cfg.conv)
  kernel = [cfg.conv(:)'./sum(cfg.conv)];
  if ~rem(length(kernel),2)
    kernel = [kernel 0];
  end
  dat = convn(dat, kernel, 'same');
end
if strcmp(cfg.derivative, 'yes'),
  dat = preproc_derivative(dat, 1, 'end');
end
if strcmp(cfg.absdiff, 'yes'),
  % this implements abs(diff(data), which is required for jump detection
  dat = abs([diff(dat, 1, 2) zeros(size(dat,1),1)]);
end
if strcmp(cfg.standardize, 'yes'),
  dat = preproc_standardize(dat, 1, size(dat,2));
end
if ~isempty(cfg.subspace),
  dat = preproc_subspace(dat, cfg.subspace);
end
if ~isempty(cfg.precision)
  % convert the data to another numeric precision, i.e. double, single or int32
  dat = cast(dat, cfg.precision);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the filter padding and do the preprocessing on the remaining trial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if begpadding~=0 || endpadding~=0
  dat = dat(:, (1+begpadding):(end-endpadding));
  if strcmp(cfg.blc, 'yes') || nargout>2
    time = time((1+begpadding):(end-endpadding));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subfunction which does the polyremoval on the data without filter-padding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dat = polyremoval(dat, polyorder, begsample, endsample);
nsamples = size(dat,2);
basis = [1:nsamples];
x = zeros(polyorder+1,nsamples);
for i = 0:polyorder
  x(i+1,:) = basis.^(i);
end
a = dat(:,begsample:endsample)/x(:,begsample:endsample);
dat = dat - a*x;
