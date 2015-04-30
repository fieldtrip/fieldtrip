function freq = ft_datatype_freq(freq, varargin)

% FT_DATATYPE_FREQ describes the FieldTrip MATLAB structure for freq data
%
% The freq data structure represents frequency or time-frequency decomposed
% channel-level data. This data structure is usually generated with the
% FT_FREQANALYSIS function.
%
% An example of a freq structure containing the powerspectrum for 306 channels
% and 102 frequencies is
%
%       dimord: 'chan_freq'          defines how the numeric data should be interpreted
%    powspctrm: [306x120 double]     the power spectum
%        label: {306x1 cell}         the channel labels
%         freq: [1x120 double]       the frequencies expressed in Hz
%          cfg: [1x1 struct]         the configuration used by the function that generated this data structure
%
% An example of a freq structure containing the time-frequency resolved
% spectral estimates of power (i.e. TFR) for 306 channels, 120 frequencies
% and 60 timepoints is
%
%       dimord: 'chan_freq_time'     defines how the numeric data should be interpreted
%    powspctrm: [306x120x60 double]  the power spectum
%        label: {306x1 cell}         the channel labels
%         freq: [1x120 double]       the frequencies, expressed in Hz
%         time: [1x60 double]        the time, expressed in seconds
%          cfg: [1x1 struct]         the configuration used by the function that generated this data structure
%
% Required fields:
%   - label, dimord, freq
%
% Optional fields:
%   - powspctrm, fouriesspctrm, csdspctrm, cohspctrm, time, labelcmb, grad, elec, cumsumcnt, cumtapcnt, trialinfo 
%
% Deprecated fields:
%   - <none>
%
% Obsoleted fields:
%   - <none>
%
% Historical fields:
%   - cfg, crsspctrm, cumsumcnt, cumtapcnt, dimord, elec, foi,
%   fourierspctrm, freq, grad, label, labelcmb, powspctrm, time, toi, see
%   bug2513
%
% Revision history:
%
% (2011/latest) The description of the sensors has changed, see FT_DATATYPE_SENS
% for further information.
%
% (2008) The presence of labelcmb in case of crsspctrm became optional,
% from now on the crsspctrm can also be represented as Nchan * Nchan.
%
% (2006) The fourierspctrm field was added as alternative to powspctrm and
% crsspctrm. The fields foi and toi were renamed to freq and time.
%
% (2003v2) The fields sgn and sgncmb were renamed into label and labelcmb.
%
% (2003v1) The initial version was defined.
%
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2011, Robert Oostenveld
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

% get the optional input arguments, which should be specified as key-value pairs
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2011';
end

if isempty(freq)
  return;
end

% ensure consistency between the dimord string and the axes that describe the data dimensions
freq = fixdimord(freq);

if ~isrow(freq.freq)
  freq.freq = freq.freq';
end
if isfield(freq, 'label') && ~iscolumn(freq.label)
  % this is not present if the dimord is chancmb_freq or chancmb_freq_time
  freq.label = freq.label';
end
if isfield(freq, 'time') && ~isrow(freq.time)
  freq.time = freq.time';
end
if ~isfield(freq, 'label') && ~isfield(freq, 'labelcmb')
  warning('data structure is incorrect since it has no channel labels');
end

switch version
  case '2011'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(freq, 'grad')
      % ensure that the gradiometer structure is up to date
      freq.grad = ft_datatype_sens(freq.grad);
    end
    
    if isfield(freq, 'elec')
      % ensure that the electrode structure is up to date
      freq.elec = ft_datatype_sens(freq.elec);
    end
 
    if isfield(freq, 'foi') && ~isfield(freq, 'freq')
      % this was still the case in early 2006
      freq.freq = freq.foi;
      freq = rmfield(freq, 'foi');
    end

    if isfield(freq, 'toi') && ~isfield(freq, 'time')
      % this was still the case in early 2006
      freq.time = freq.toi;
      freq = rmfield(freq, 'toi');
    end

  case '2008'
    % there are no known conversions for backward or forward compatibility support

  case '2006'
    % there are no known conversions for backward or forward compatibility support

  case '2003v2'
    % there are no known conversions for backward or forward compatibility support

  case '2003v1'
    % there are no known conversions for backward or forward compatibility support

  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for freq datatype', version);
end

