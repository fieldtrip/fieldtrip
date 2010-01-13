function [freq] = freqanalysis(cfg, data);

% FREQANALYSIS performs frequency and time-frequency analysis
% on time series data over multiple trials
%
% Use as
%   [freq] = freqanalysis(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the PREPROCESSING function. The configuration depends on the type
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
% read the help of the respective subfunction FREQANALYSIS_XXX for the
% corresponding parameter options and for a detailed explanation of each method.
%
% See also FREQANALYSIS_MTMFFT, FREQANALYSIS_MTMCONVOL, FREQANALYSIS_MTMWELCH
% FREQANALYSIS_WLTCONVOL, FREQANALYSIS_TFR

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

% call the corresponding function
[freq] = feval(sprintf('ft_freqanalysis_%s',lower(cfg.method)), cfg, data);
