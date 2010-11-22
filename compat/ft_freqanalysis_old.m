function [freq] = ft_freqanalysis_old(cfg, data)

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
% cfg.correctt_ftimwin (set to yes to try to determine new t_ftimwins based on correct cfg.foi)
% cfg.channel
% cfg.channelcmb
% cfg.inputfile  = one can specifiy preanalysed saved data as input
% cfg.outputfile = one can specify output as file to save to disk

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
% $Id: ft_freqanalysis.m 2145 2010-11-22 13:28:48Z jansch $

fieldtripdefs

% deprecated
warning('this version of ft_freqanalysis is deprecated, it is no longer being maintained')

% defaults for optional input/ouputfile
if ~isfield(cfg, 'inputfile'),  cfg.inputfile                   = [];    end
if ~isfield(cfg, 'outputfile'), cfg.outputfile                  = [];    end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
  end
end


% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw', 'comp', 'mvar'}, 'feedback', 'yes', 'hasoffset', 'yes', 'hastrialdef', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'renamed',     {'label', 'channel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'sgn',   'channel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'labelcmb', 'channelcmb'});
cfg = ft_checkconfig(cfg, 'renamed',     {'sgncmb',   'channelcmb'});
cfg = ft_checkconfig(cfg, 'required',    {'method'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'fft',    'mtmfft'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'convol', 'mtmconvol'});

% select trials of interest
if ~isfield(cfg, 'trials'),   cfg.trials = 'all';  end % set the default
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = ft_selectdata(data, 'rpt', cfg.trials);
end

% Switch for the old implementation (will be in ft_freqanalysis.m in a while)
switch cfg.method
  case 'mtmconvol'
    cfg.method = [cfg.method '_old'];
  case 'mtmfft'
    cfg.method = [cfg.method '_old'];
  case 'mtmwelch'
    cfg.method = [cfg.method '_old'];
  case 'wltconvol'
    cfg = ft_checkconfig(cfg, 'renamedval', {'method','wltconvol','wavelet'});
    cfg.method = [cfg.method '_old'];
  case 'tfr'
    cfg = ft_checkconfig(cfg, 'renamedval', {'method','tfr','convol'});
    cfg.method = [cfg.method '_old'];
  case 'mvar'
    cfg.method = [cfg.method '_old'];
end










% call old functions
[freq] = feval(sprintf('ft_freqanalysis_%s',lower(cfg.method)), cfg, data);
  




% copy the trial specific information into the output
if isfield(data, 'trialinfo'),
  freq.trialinfo = data.trialinfo;
  % FIXME this strictly only allowed for frequency data with repetitions
end

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'freq', freq); % use the variable name "data" in the output file
end










