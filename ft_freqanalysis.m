function [freq] =ft_freqanalysis(cfg, data);

% FT_FREQANALYSIS performs frequency and time-frequency analysis
% on time series data over multiple trials
%
% Use as
%   [freq] = ft_freqanalysis(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration depends on the type
% of computation that you want to perform.
%
% The configuration should contain:
%   cfg.method     = different methods of calculating the spectra
%                    'mtmfft', analyses an entire spectrum for the entire data
%                     length, implements multitaper frequency transformation
%                    'mtmconvol', implements multitaper time-frequency transformation
%                     based on multiplication in the frequency domain
%                    'mtmwelch', performs frequency analysis using Welch's averaged
%                     modified periodogram method of spectral estimation
%                    'wltconvol', implements wavelet time frequency transformation
%                     (using Morlet wavelets) based on multiplication in the frequency domain
%                    'tfr', implements wavelet time frequency transformation
%                     (using Morlet wavelets) based on convolution in the time domain
%
%
% The other cfg options depend on the method that you select. You should
% read the help of the respective subfunction FT_FREQANALYSIS_XXX for the
% corresponding parameter options and for a detailed explanation of each method.
%
% See also FT_FREQANALYSIS_MTMFFT, FT_FREQANALYSIS_MTMCONVOL, FT_FREQANALYSIS_MTMWELCH
% FT_FREQANALYSIS_WLTCONVOL, FT_FREQANALYSIS_TFR

% Undocumented local options:
% cfg.label
% cfg.labelcmb
% cfg.sgn
% cfg.sgncmb

% Copyright (C) 2003-2006, F.C. Donders Centre, Pascal Fries
% Copyright (C) 2004-2006, F.C. Donders Centre, Markus Siegel
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

% check if the input data is valid for this function
data = checkdata(data, 'datatype', {'raw', 'comp', 'mvar'}, 'feedback', 'yes', 'hasoffset', 'yes');

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'renamed',     {'label', 'channel'});
cfg = checkconfig(cfg, 'renamed',     {'sgn',   'channel'});
cfg = checkconfig(cfg, 'renamed',     {'labelcmb', 'channelcmb'});
cfg = checkconfig(cfg, 'renamed',     {'sgncmb',   'channelcmb'});
cfg = checkconfig(cfg, 'required',    {'method'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'fft',    'mtmfft'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'convol', 'mtmconvol'});

% select trials of interest
if ~isfield(cfg, 'trials'),   cfg.trials = 'all';  end % set the default
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = selectdata(data, 'rpt', cfg.trials);
  if isfield(data, 'cfg') % try to locate the trl in the nested configuration
    cfg.trlold = findcfg(data.cfg, 'trlold');
    cfg.trl    = findcfg(data.cfg, 'trl');
  end
end

if true
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % HERE THE OLD IMPLEMENTATION STARTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [freq] = feval(sprintf('ft_freqanalysis_%s',lower(cfg.method)), cfg, data);

else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % HERE THE NEW IMPLEMENTATION STARTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % do the bookkeeping that is required for the input data
  % chansel = ...;

  % do the bookkeeping that is required for the computation
  offset    = data.offset;
  trllength = zeros(1, ntrials);

  for trllop=1:ntrials
    trllength(i) = size(data.trial{trllop}, 2);
  end

  if strcmp(cfg.padding, 'maxperlen')
    padding = max(trllength);
    cfg.padding = padding/data.fsample;
  else
    padding = cfg.padding*data.fsample;
    if padding<max(trllength)
      error('the specified padding is too short');
    end
  end

  % these don't change over trials
  options = {'fsample', data.fsample, 'padding', padding};

  % do the bookkeeping that is required for the output data
  % ...

  state = [];
  for trllop=1:ntrials

    dat = data.trial{trllop}(chansel,:);

    % do the spectral decompisition of this trial
    switch cfg.method
      case 'mtmfft'
        [spectrum, state] = specest_mtmfft(   dat, options{:}, 'state', state, 'offset', offset(trllop));
      case 'mtmconvol'
        [spectrum, state] = specest_mtmconvol(dat, options{:}, 'state', state, 'offset', offset(trllop));
      case 'wltconvol'
        [spectrum, state] = specest_wltconvol(dat, options{:}, 'state', state, 'offset', offset(trllop));
      otherwise
        error('method %s is unknown', cfg.method);
    end % switch

    % do the bookkeping to fit the spectral decomposition in the correct
    % output representation

    % ...

  end % for ntrials

end % if old or new implementation

