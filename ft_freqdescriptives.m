function [output] = ft_freqdescriptives(cfg, freq)

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
%   cfg.foilim        = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%   cfg.toilim        = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
%
% A variance estimate can only be computed if results from trials and/or
% tapers have been kept.
%
% Descriptive statistics of bivariate metrics is not computed by this function anymore. To this end you
% should use FT_CONNECTIVITYANALYSIS.
%
% See also FT_FREQANALYSIS, FT_FREQSTATISTICS, FT_FREQBASELINE, FT_CONNECTIVITYANALYSIS

% Undocumented local options:
% cfg.feedback
% cfg.latency
% cfg.previous
% cfg.version
% cfg.inputfile  = one can specifiy preanalysed saved data as input
% cfg.outputfile = one can specify output as file to save to disk

% Copyright (C) 2004-2006, Pascal Fries & Jan-Mathijs Schoffelen, F.C. Donders Centre
% Copyright (C) 2010, Jan-Mathijs Schoffelen, F.C. Donders Centre
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'renamed',     {'jacknife',   'jackknife'});

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

% set the defaults
if ~isfield(cfg, 'feedback'),    cfg.feedback      = 'textbar'; end
if ~isfield(cfg, 'jackknife'),   cfg.jackknife     = 'no';      end
if ~isfield(cfg, 'variance'),    cfg.variance      = 'no';      end
if ~isfield(cfg, 'trials'),      cfg.trials        = 'all';     end
if ~isfield(cfg, 'channel'),     cfg.channel       = 'all';     end
if ~isfield(cfg, 'foilim'),      cfg.foilim        = 'all';     end
if ~isfield(cfg, 'toilim'),      cfg.toilim        = 'all';     end
if ~isfield(cfg, 'keeptrials'),  cfg.keeptrials    = 'no';      end
if ~isfield(cfg, 'inputfile'),   cfg.inputfile     = [];        end
if ~isfield(cfg, 'outputfile'),  cfg.outputfile    = [];        end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    freq = loadvar(cfg.inputfile, 'freq');
  end
end

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'datatype', {'freq', 'freqmvar'}, 'feedback', 'yes');

% determine some specific details of the input data
hasrpt   = ~isempty(strfind(freq.dimord, 'rpt')) || ~isempty(strfind(freq.dimord, 'subj'));
hastim   = ~isempty(strfind(freq.dimord, 'time'));

varflg   = strcmp(cfg.variance,  'yes');
jckflg   = strcmp(cfg.jackknife, 'yes');
keepflg  = strcmp(cfg.keeptrials, 'yes');

% check sensibility of configuration
if sum([varflg keepflg]>1),               error('you should specify only one of cfg.keeptrials or cfg.variance');                                             end
if ~hasrpt && (varflg || keepflg),        error('a variance-estimate or a single trial estimate without repeated observations in the input is not possible'); end
if ~hasrpt && ~strcmp(cfg.trials, 'all'), error('trial selection requires input data with repeated observations');                                            end
if ~varflg && jckflg,                     warning('you specified cfg.jackknife = ''yes'' and cfg.variance = ''no'': no variance will be computed');           end

% select data of interest
if            ~strcmp(cfg.foilim,  'all'), freq = ft_selectdata(freq, 'foilim', cfg.foilim); end
if hastim, if ~strcmp(cfg.toilim,  'all'), freq = ft_selectdata(freq, 'toilim', cfg.toilim); end; end
if hasrpt, if ~strcmp(cfg.trials,  'all'), freq = ft_selectdata(freq, 'rpt',    cfg.trials); end; end

if ~strcmp(cfg.channel, 'all'),
  channel = ft_channelselection(cfg.channel, freq.label);
  freq    = ft_selectdata(freq, 'channel', channel);
end

% get data in the correct representation
freq = ft_checkdata(freq, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});

if jckflg,
  freq = ft_selectdata(freq, 'jackknife', 1);
end

if varflg,
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
  
  if jckflg,
    bias = (n-1).^2;
  else,
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
  powspctrm = reshape(nanmean(freq.powspctrm,1), [siz(2:end)]);
else
  %nothing to do
  powspctrm = freq.powspctrm;
end

if hasrpt && ~keepflg,
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
try, output.time      = freq.time;      end;
try, output.grad      = freq.grad;      end;
try, output.cumtapcnt = freq.cumtapcnt; end;
try, output.cumsumcnt = freq.cumsumcnt; end;
output.powspctrm      = powspctrm;
try, output.powspctrmsem = powspctrmsem; end;

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';

try, cfg.previous = freq.cfg; end

% remember the configuration details
output.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', output); % use the variable name "data" in the output file
end
