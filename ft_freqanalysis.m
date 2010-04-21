function [freq] =ft_freqanalysis(cfg, data, flag);

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

fieldtripdefs

%allow for both the new and old implementation to be changed with a flag
%input

if nargin < 3
    flag = 0;
end

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

if ~flag
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
  %offset    = data.offset;
  
  if ~isfield(cfg, 'padding'), cfg.padding = [];   end
  if ~isfield(cfg, 'output'),  cfg.output = 'pow'; end
  if ~isfield(cfg, 'taper'),   cfg.taper =  'dpss';     end
  if ~isfield(cfg, 'method'), error('you must specify a method in cfg.method'); end
  if ~isfield(cfg, 'foi'),     cfg.foi = [];       end
  if isequal(cfg.taper, 'dpss') && not(isfield(cfg, 'tapsmofrq'))
      error('you must specify a smoothing parameter with taper = dpss'); 
  end  
  
  ntrials = size(data.trial,2);
  trllength = zeros(1, ntrials);

  for trllop=1:ntrials
    trllength(trllop) = size(data.trial{trllop}, 2);
  end
  
  if not(isempty(cfg.padding))
      if strcmp(cfg.padding, 'maxperlen')
        padding = max(trllength);
        cfg.padding = padding/data.fsample;
      else
        padding = cfg.padding*data.fsample;
        if padding<max(trllength)
          error('the specified padding is too short');
        end
      end
  end
  % these don't change over trials
  options = {'pad', cfg.padding, 'taper', cfg.taper, 'freqoi', cfg.foi, 'tapsmofrq', cfg.tapsmofrq}; 
  
  % do the bookkeeping that is required for the output data
  freq  = [];
  numsmp = length(data.time{1});
  flag = 0;
  for trllop=1:ntrials

    dat = data.trial{trllop}; %chansel has already been performed 
    time = data.time{trllop};
    

    % do the spectral decompisition of this trial
    switch cfg.method
      case 'mtmfft'
        [spectrum, foi] = specest_mtmfft(dat, time, options{:}); %Add offset option later 'offset', offset(trllop));
        if flag == 0
            fourierspctrm = zeros(ntrials,size(spectrum,1),size(spectrum,2),size(spectrum,3)); %this matrix is trials by tapers by channel by frequency
            flag = 1;
        end
      case 'mtmconvol'
        [spectrum, foi] = specest_mtmconvol(dat, time, options{:});
      case 'wltconvol'
        [spectrum, foi] = specest_wltconvol(dat, time, options{:});
      otherwise
        error('method %s is unknown', cfg.method);
    end % switch

    fourierspctrm(trllop,:,:,:) = spectrum;

  end % for ntrials
  %now get the output in the correct format
  if isequal(cfg.output, 'pow')
    freq.powspctrm = 2 .* (fourierspctrm .* conj(fourierspctrm)) ./ numsmp; %cf Numercial Receipes 13.4.9
    freq.powspctrm = reshape(freq.powspctrm, size(freq.powspctrm,1)*size(freq.powspctrm,2),size(freq.powspctrm,3),size(freq.powspctrm,4));
    freq.dimord = 'rpttap_chan_freq';
    freq.freq = foi;
  elseif isequal(cfg.output, 'fourier')
    freq.fourierspctrm = (fourierspctrm) .* sqrt(2 ./ numsmp); %cf Numercial Receipes 13.4.9
    freq.fourierspctrm = reshape(freq.fourierspctrm, size(freq.fourierspctrm,1)*size(freq.fourierspctrm,2),size(freq.fourierspctrm,3),size(freq.fourierspctrm,4));
    freq.dimord = 'rpttap_chan_freq'; 
    freq.freq = foi;
  elseif isequal(cfg.output, 'csd')
    freq.dimord = 'rpttap_chan_freq';
    freq.freq = foi;
    freq.label = data.label;
    freq.cumtapcnt = size(spectrum,1)*zeros(ntrials,1);
    freq.fourierspctrm = (fourierspctrm) .* sqrt(2 ./ numsmp); %cf Numercial Receipes 13.4.9
    freq.fourierspctrm = reshape(freq.fourierspctrm, size(freq.fourierspctrm,1)*size(freq.fourierspctrm,2),size(freq.fourierspctrm,3),size(freq.fourierspctrm,4));
    freq = fixcsd(freq, 'full', []);
  else
    error('output is not recognized',cfg.output);
  end  

end % if old or new implementation

