function [output] = ft_freqdescriptives(cfg, freq)

% FT_FREQDESCRIPTIVES computes descriptive univariate statistics of 
% the frequency or time-frequency decomposition of the EEG/MEG signal, 
% thus the powerspectrum and its standard error.
%
% Use as
%   [freq] = ft_freqdescriptives(cfg, freq)
%
% The data in freq should be organised in a structure as obtained from
% from the FT_FREQANALYSIS function. The output structure is comparable
% to the input structure and can be used in most functions that require
% a freq input.
%
% The configuration options are
%   cfg.variance      = 'yes' or 'no', estimate standard error in the standard way (default = 'no')
%   cfg.jackknife     = 'yes' or 'no', estimate standard error by means of the jack-knife (default = 'no')
%   cfg.keeptrials    = 'yes' or 'no', estimate single trial power (useful for fourier data) (default = 'no')
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                       see CHANNELSELECTION for details
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

% Copyright (C) 2004-2006, Pascal Fries & Jan-Mathijs Schoffelen, F.C. Donders Centre
% Copyright (C) 2010, Jan-Mathijs Schoffelen, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailed information

fieldtripdefs

% check if the input data is valid for this function
freq = checkdata(freq, 'datatype', {'freq', 'freqmvar'}, 'feedback', 'yes');

% throw warnings for the deprecated options
cfg = checkconfig(cfg, 'deprecated', 'biascorrect');
cfg = checkconfig(cfg, 'deprecated', 'channelcmb');
cfg = checkconfig(cfg, 'deprecated', 'cohmethod');
cfg = checkconfig(cfg, 'deprecated', 'combinemethod');
cfg = checkconfig(cfg, 'deprecated', 'complex');
cfg = checkconfig(cfg, 'deprecated', 'combinechan');
cfg = checkconfig(cfg, 'deprecated', 'keepfourier');
cfg = checkconfig(cfg, 'deprecated', 'partchan');
cfg = checkconfig(cfg, 'deprecated', 'pseudovalue');

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'renamed',     {'jacknife',   'jackknife'});

% determine some specific details of the input data
hasrpt   = ~isempty(strfind(freq.dimord, 'rpt')) || ~isempty(strfind(freq.dimord, 'subj'));
hastim   = ~isempty(strfind(freq.dimord, 'time'));

% set the defaults
if ~isfield(cfg, 'feedback'),    cfg.feedback      = 'textbar'; end
if ~isfield(cfg, 'jackknife'),   cfg.jackknife     = 'no';      end
if ~isfield(cfg, 'variance'),    cfg.variance      = 'no';      end
if ~isfield(cfg, 'trials'),      cfg.trials        = 'all';     end
if ~isfield(cfg, 'channel'),     cfg.channel       = 'all';     end
if ~isfield(cfg, 'foilim'),      cfg.foilim        = 'all';     end
if ~isfield(cfg, 'toilim'),      cfg.toilim        = 'all';     end
if ~isfield(cfg, 'keeptrials'),  cfg.keeptrials    = 'no';      end

varflg   = strcmp(cfg.variance,  'yes');
jckflg   = strcmp(cfg.jackknife, 'yes');
keepflg  = strcmp(cfg.keeptrials, 'yes');

% check sensibility of configuration
if sum([varflg keepflg]>1),               error('you should specify only one of cfg.keeptrials or cfg.variance');                                             end
if ~hasrpt && (varflg || keepflg),        error('a variance-estimate or a single trial estimate without repeated observations in the input is not possible'); end
if ~hasrpt && ~strcmp(cfg.trials, 'all'), error('trial selection requires input data with repeated observations');                                            end
if ~varflg && jckflg,                     warning('you specified cfg.jackknife = ''yes'' and cfg.variance = ''no'': no variance will be computed');           end

% select data of interest
if            ~strcmp(cfg.foilim,  'all'), freq = selectdata(freq, 'foilim', cfg.foilim); end
if hastim, if ~strcmp(cfg.toilim,  'all'), freq = selectdata(freq, 'toilim', cfg.toilim); end; end
if hasrpt, if ~strcmp(cfg.trials,  'all'), freq = selectdata(freq, 'rpt',    cfg.trials); end; end

if ~strcmp(cfg.channel, 'all'),
  channel = channelselection(cfg.channel, freq.label);
  freq    = selectdata(freq, 'channel', channel);
end

% get data in the correct representation
freq = checkdata(freq, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});

if jckflg,
  freq = selectdata(freq, 'jackknife', 1);
end

if varflg,
  siz    = [size(freq.powspctrm) 1];
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  n      = zeros(siz(2:end));
  progress('init', cfg.feedback, 'computing power...');
  for j = 1:siz(1)
    progress(j/siz(1), 'computing power for replicate %d from %d\n', j, siz(1));
    tmp    = reshape(freq.powspctrm(j,:,:,:), siz(2:end));
    n      = n + double(isfinite(tmp));
    tmp(~isfinite(tmp)) = 0;
    outsum = outsum + tmp;
    outssq = outssq + tmp.^2;
  end 
  progress('close');

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

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

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
