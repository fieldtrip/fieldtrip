function [freq] = ft_freqdescriptives(cfg, freq)

% FT_FREQDESCRIPTIVES computes descriptive univariate statistics of
% the frequency or time-frequency decomposition of the EEG/MEG signal,
% thus the powerspectrum and its standard error.
%
% Use as
%   [freq] = ft_freqdescriptives(cfg, freq)
%   [freq] = ft_freqdescriptives(cfg, freqmvar)
%
% The data in freq should be organised in a structure as obtained from
% from the FT_FREQANALYSIS or FT_MVARANALYSIS function. The output structure is comparable
% to the input structure and can be used in most functions that require
% a freq input.
%
% The configuration options are
%   cfg.variance      = 'yes' or 'no', estimate standard error in the standard way (default = 'no')
%   cfg.jackknife     = 'yes' or 'no', estimate standard error by means of the jack-knife (default = 'no')
%   cfg.keeptrials    = 'yes' or 'no', estimate single trial power (useful for fourier data) (default = 'no')
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                       see FT_CHANNELSELECTION for details
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.frequency     = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%   cfg.latency       = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
%
% A variance estimate can only be computed if results from trials and/or
% tapers have been kept.
%
% Descriptive statistics of bivariate metrics is not computed by this function anymore. To this end you
% should use FT_CONNECTIVITYANALYSIS.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_FREQANALYSIS, FT_FREQSTATISTICS, FT_FREQBASELINE, FT_CONNECTIVITYANALYSIS

% Undocumented local options:
% cfg.feedback
% cfg.latency
% cfg.previous
% cfg.version

% Copyright (C) 2004-2006, Pascal Fries & Jan-Mathijs Schoffelen
% Copyright (C) 2010, Jan-Mathijs Schoffelen
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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'datatype', {'freq', 'freqmvar'}, 'feedback', 'yes');
% get data in the correct representation, it should only have power
freq = ft_checkdata(freq, 'cmbstyle', 'sparsewithpow', 'channelcmb', {});

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels', 'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamed',    {'jacknife', 'jackknife'});

% throw warnings for the deprecated options
cfg = ft_checkconfig(cfg, 'deprecated', 'biascorrect');
cfg = ft_checkconfig(cfg, 'deprecated', 'channelcmb');
cfg = ft_checkconfig(cfg, 'deprecated', 'cohmethod');
cfg = ft_checkconfig(cfg, 'deprecated', 'combinemethod');
cfg = ft_checkconfig(cfg, 'deprecated', 'complex');
cfg = ft_checkconfig(cfg, 'deprecated', 'combinechan');
cfg = ft_checkconfig(cfg, 'deprecated', 'keepfourier');
cfg = ft_checkconfig(cfg, 'deprecated', 'partchan');
cfg = ft_checkconfig(cfg, 'deprecated', 'pseudovalue');
cfg = ft_checkconfig(cfg, 'renamed', {'toilim' 'latency'});
cfg = ft_checkconfig(cfg, 'renamed', {'foilim' 'frequency'});

% set the defaults
cfg.feedback   = ft_getopt(cfg, 'feedback',   'textbar');
cfg.jackknife  = ft_getopt(cfg, 'jackknife',  'no');
cfg.variance   = ft_getopt(cfg, 'variance',   'no');
cfg.trials     = ft_getopt(cfg, 'trials',     'all', 1);
cfg.channel    = ft_getopt(cfg, 'channel',    'all');
cfg.frequency  = ft_getopt(cfg, 'frequency',  'all');
cfg.latency    = ft_getopt(cfg, 'latency',    'all');
cfg.keeptrials = ft_getopt(cfg, 'keeptrials', 'no');

% determine some specific details of the input data
hasrpt   = ~isempty(strfind(freq.dimord, 'rpt')) || ~isempty(strfind(freq.dimord, 'subj'));
hastim   = ~isempty(strfind(freq.dimord, 'time'));

varflg   = strcmp(cfg.variance,  'yes');
jckflg   = strcmp(cfg.jackknife, 'yes');
keepflg  = strcmp(cfg.keeptrials, 'yes');

% check sensibility of configuration
if sum([varflg keepflg]>1),               ft_error('you should specify only one of cfg.keeptrials or cfg.variance');                                             end
if ~hasrpt && (varflg || keepflg),        ft_error('a variance-estimate or a single trial estimate without repeated observations in the input is not possible'); end
if ~hasrpt && ~strcmp(cfg.trials, 'all'), ft_error('trial selection requires input data with repeated observations');                                            end
if ~varflg && jckflg,                     varflg = 1; end

% select data of interest
tmpcfg = keepfields(cfg, {'trials', 'channel', 'latency', 'frequency', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
freq = ft_selectdata(tmpcfg, freq);
% restore the provenance information
[cfg, freq] = rollback_provenance(cfg, freq);

if jckflg
  % the data is 'sparsewithpow', so it contains a powspctrm and optionally a crsspctrm
  % the checking of a 'rpt' is handled above, so it can be assumed that the 'rpt' is the
  % first dimension
  nrpt           = size(freq.powspctrm,1);
  sumpowspctrm   = sum(freq.powspctrm,1);
  freq.powspctrm = (sumpowspctrm(ones(nrpt,1),:,:,:,:) - freq.powspctrm)./(nrpt-1);
  clear sumpowspctrm;
  if isfield(freq, 'crsspctrm')
    sumcrsspctrm   = sum(freq.crsspctrm,1);
    freq.crsspctrm = (sumcrsspctrm(ones(nrpt,1),:,:,:,:) - freq.crsspctrm)./(nrpt-1);
    clear sumcrsspctrm;
  end
end

if varflg
  siz    = [size(freq.powspctrm) 1];
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  n      = zeros(siz(2:end));
  ft_progress('init', cfg.feedback, 'computing power...');
  for j = 1:siz(1)
    ft_progress(j/siz(1), 'computing power for replicate %d from %d\n', j, siz(1));
    tmp    = reshape(freq.powspctrm(j,:,:,:), siz(2:end));
    n      = n + double(isfinite(tmp));
    tmp(~isfinite(tmp)) = 0;
    outsum = outsum + tmp;
    outssq = outssq + tmp.^2;
  end
  ft_progress('close');

  if jckflg
    bias = (n-1).^2;
  else
    bias = 1;
  end

  powspctrm    = outsum./n;
  powspctrmsem = sqrt(bias.*(outssq - (outsum.^2)./n)./(n - 1)./n);
elseif keepflg
  %nothing to do
  powspctrm = freq.powspctrm;
elseif hasrpt
  %compute average only
  siz       = [size(freq.powspctrm) 1];
  powspctrm = reshape(nanmean(freq.powspctrm,1), siz(2:end));
else
  %nothing to do
  powspctrm = freq.powspctrm;
end

if hasrpt && ~keepflg
  dimtok    = tokenize(freq.dimord, '_');
  newdimord = dimtok{2};
  for k = 3:numel(dimtok)
    newdimord = [newdimord,'_',dimtok{k}];
  end
else
  newdimord = freq.dimord;
end

% create the output-structure
output                = [];
output.dimord         = newdimord;
output.freq           = freq.freq;
output.label          = freq.label;
output.powspctrm      = powspctrm;
if isfield(freq, 'time'),      output.time      = freq.time;      end
if isfield(freq, 'grad'),      output.grad      = freq.grad;      end
if isfield(freq, 'cumtapcnt'), output.cumtapcnt = freq.cumtapcnt; end
if isfield(freq, 'cumsumcnt'), output.cumsumcnt = freq.cumsumcnt; end
if exist('powspctrmsem', 'var'), output.powspctrmsem = powspctrmsem; end

% remember the trialinfo
if strcmp(cfg.keeptrials, 'yes') && isfield(freq, 'trialinfo')
  output.trialinfo = freq.trialinfo;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous freq

% rename the output variable to accomodate the savevar postamble
freq = output;

ft_postamble provenance freq
ft_postamble history    freq
ft_postamble savevar    freq
