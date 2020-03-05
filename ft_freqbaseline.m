function [freq] = ft_freqbaseline(cfg, freq)

% FT_FREQBASELINE performs baseline normalization for time-frequency data
%
% Use as
%    [freq] = ft_freqbaseline(cfg, freq)
% where the freq data comes from FT_FREQANALYSIS and the configuration
% should contain
%   cfg.baseline     = [begin end] (default = 'no'), alternatively an
%                      Nfreq x 2 matrix can be specified, that provides
%                      frequency specific baseline windows.
%   cfg.baselinetype = 'absolute', 'relative', 'relchange', 'normchange', 'db', 'vssum' or 'zscore' (default = 'absolute')
%   cfg.parameter    = field for which to apply baseline normalization, or
%                      cell-array of strings to specify multiple fields to normalize
%                      (default = 'powspctrm')
%
% See also FT_FREQANALYSIS, FT_TIMELOCKBASELINE, FT_FREQCOMPARISON,
% FT_FREQGRANDAVERAGE

% Undocumented local options:
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%   cfg.outputfile = one can specify output as file to save to disk

% Copyright (C) 2004-2006, Marcel Bastiaansen
% Copyright (C) 2005-2006, Robert Oostenveld
% Copyright (C) 2011, Eelke Spaak
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
ft_preamble debug
ft_preamble loadvar freq
ft_preamble provenance freq
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'datatype', {'freq+comp', 'freq'}, 'feedback', 'yes');

% update configuration fieldnames
cfg              = ft_checkconfig(cfg, 'renamed', {'param', 'parameter'});

% set the defaults
cfg.baseline     =  ft_getopt(cfg, 'baseline', 'no');
cfg.baselinetype =  ft_getopt(cfg, 'baselinetype', 'absolute');
cfg.parameter    =  ft_getopt(cfg, 'parameter', 'powspctrm');

% check validity of input options
cfg =               ft_checkopt(cfg, 'baseline', {'char', 'doublevector', 'doublematrix'});
cfg =               ft_checkopt(cfg, 'baselinetype', 'char', {'absolute', 'relative', 'relchange', 'normchange', 'db', 'vssum','zscore'});
cfg =               ft_checkopt(cfg, 'parameter', {'char', 'charcell'});

% make sure cfg.parameter is a cell-array of strings
if (~isa(cfg.parameter, 'cell'))
  cfg.parameter = {cfg.parameter};
end

% is input consistent?
if ischar(cfg.baseline) && strcmp(cfg.baseline, 'no') && ~isempty(cfg.baselinetype)
  ft_warning('no baseline correction done');
end

% process possible yes/no value of cfg.baseline
if ischar(cfg.baseline) && strcmp(cfg.baseline, 'yes')
  % default is to take the prestimulus interval
  cfg.baseline = [-inf 0];
elseif ischar(cfg.baseline) && strcmp(cfg.baseline, 'no')
  % nothing to do
  return
end

% allow for baseline to be nfreq x 2
if size(cfg.baseline,1)==numel(freq.freq) && size(cfg.baseline,2)==2
  % this is ok
elseif numel(cfg.baseline)==2
  % this is also ok
  cfg.baseline = cfg.baseline(:)'; % ensure row vector
else
  ft_error('cfg.baseline should either be a string, a 1x2 vector, or an Nfreqx2 matrix');
end

% check if the field of interest is present in the data
if (~all(isfield(freq, cfg.parameter)))
  ft_error('cfg.parameter should be a string or cell-array of strings referring to (a) field(s) in the freq input structure')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output structure
freqOut = keepfields(freq, {'label' 'freq' 'dimord' 'time'});
freqOut = copyfields(freq, freqOut,...
  {'grad', 'elec', 'trialinfo', 'topo', 'topolabel', 'unmixing'});

% loop over all fields that should be normalized
for k = 1:numel(cfg.parameter)
  par = cfg.parameter{k};

  if strcmp(freq.dimord, 'chan_freq_time')

    freqOut.(par) = ...
      performNormalization(freq.time, freq.(par), cfg.baseline, cfg.baselinetype);

  elseif strcmp(freq.dimord, 'rpt_chan_freq_time') || strcmp(freq.dimord, 'chan_chan_freq_time') || strcmp(freq.dimord, 'subj_chan_freq_time')

    freqOut.(par) = zeros(size(freq.(par)));

    % loop over trials, perform normalization per trial
    for l = 1:size(freq.(par), 1)
      tfdata = freq.(par)(l,:,:,:);
      siz    = size(tfdata);
      tfdata = reshape(tfdata, siz(2:end));
      freqOut.(par)(l,:,:,:) = ...
        performNormalization(freq.time, tfdata, cfg.baseline, cfg.baselinetype);
    end

  else
    ft_error('unsupported data dimensions: %s', freq.dimord);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output scaffolding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(cfg.parameter)==1
  % convert from cell-array to string
  cfg.parameter = cfg.parameter{1};
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   freq

% rename the output variable to accomodate the savevar postamble
freq = freqOut;

ft_postamble provenance freq
ft_postamble history    freq
ft_postamble savevar    freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that actually performs the normalization on an arbitrary quantity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = performNormalization(timeVec, data, baseline, baselinetype)

baselineTimes = false(size(baseline,1),numel(timeVec));
for k = 1:size(baseline,1)
  baselineTimes(k,:) = (timeVec >= baseline(k,1) & timeVec <= baseline(k,2));
end

if length(size(data)) ~= 3
  ft_error('time-frequency matrix should have three dimensions (chan,freq,time)');
end

% compute mean of time/frequency quantity in the baseline interval,
% ignoring NaNs, and replicate this over time dimension
if size(baselineTimes,1)==size(data,2)
  % do frequency specific baseline
  meanVals = nan+zeros(size(data));
  for k = 1:size(baselineTimes,1)
    meanVals(:,k,:) = repmat(nanmean(data(:,k,baselineTimes(k,:)), 3), [1 1 size(data, 3)]);
  end
else
  meanVals = repmat(nanmean(data(:,:,baselineTimes), 3), [1 1 size(data, 3)]);
end

if (strcmp(baselinetype, 'absolute'))
  data = data - meanVals;
elseif (strcmp(baselinetype, 'relative'))
  data = data ./ meanVals;
elseif (strcmp(baselinetype, 'relchange'))
  data = (data - meanVals) ./ meanVals;
elseif (strcmp(baselinetype, 'normchange')) || (strcmp(baselinetype, 'vssum'))
  data = (data - meanVals) ./ (data + meanVals);
elseif (strcmp(baselinetype, 'db'))
  data = 10*log10(data ./ meanVals);
elseif (strcmp(baselinetype,'zscore'))
    stdVals = repmat(nanstd(data(:,:,baselineTimes),1, 3), [1 1 size(data, 3)]);
    data=(data-meanVals)./stdVals;
else
  ft_error('unsupported method for baseline normalization: %s', baselinetype);
end
