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
%   - powspctrm, fouriesspctrm, csdspctrm, cohspctrm, time, labelcmb, grad, elec
%
% Deprecated fields:
%   - <none>
%
% Obsoleted fields:
%   - <none>
%
% Revision history:
%
% (2008/latest) The presence of labelcmb in case of crsspctrm became optional,
% from now on the crsspctrm can also be represented as Nchan * Nchan.
%
% (2006) The fourierspctrm field was added as alternative to powspctrm and
% crsspctrm.
%
% (2003) The initial version was defined.
%
% See also FT_DATATYPE and FT_DATATYPE_xxx

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
version = keyval('version', varargin); if isempty(version), version = 'latest'; end

if strcmp(version, 'latest')
  version = '2008';
end

% ensure consistency between the dimord string and the axes that describe the data dimensions
freq = fixdimord(freq);

switch version
  case '2008'
    % there are no known conversions for backward or forward compatibility support

  case '2006'
    % there are no known conversions for backward or forward compatibility support

  case '2003'
    % there are no known conversions for backward or forward compatibility support

  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for freq datatype', version);
end

