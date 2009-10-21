function [stat] = freqstatistics(cfg, varargin)

% FREQSTATISTICS computes significance probabilities and/or critical values of a parametric statistical test 
% or a non-parametric permutation test.
%
% Use as
%   [stat] = freqstatistics(cfg, freq1, freq2, ...)
% where the input data is the result from FREQANALYSIS, FREQDESCRIPTIVES
% or from FREQGRANDAVERAGE.
%
% The configuration can contain the following options for data selection
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see CHANNELSELECTION for details
%   cfg.latency     = [begin end] in seconds or 'all' (default = 'all')
%   cfg.frequency   = [begin end], can be 'all'       (default = 'all')
%   cfg.avgoverchan = 'yes' or 'no'                   (default = 'no')
%   cfg.avgovertime = 'yes' or 'no'                   (default = 'no')
%   cfg.avgoverfreq = 'yes' or 'no'                   (default = 'no')
%   cfg.parameter   = string                          (default = 'powspctrm')
%
% Furthermore, the configuration should contain
%   cfg.method       = different methods for calculating the significance probability and/or critical value
%                    'montecarlo' get Monte-Carlo estimates of the significance probabilities and/or critical values from the permutation distribution,
%                    'analytic'   get significance probabilities and/or critical values from the analytic reference distribution (typically, the sampling distribution under the null hypothesis),
%                    'stats'      use a parametric test from the Matlab statistics toolbox,
%                    'glm'        use a general linear model approach.
%
% The other cfg options depend on the method that you select. You
% should read the help of the respective subfunction STATISTICS_XXX
% for the corresponding configuration options and for a detailed
% explanation of each method.
%
% See also FREQANALYSIS, FREQDESCRIPTIVES, FREQGRANDAVERAGE

% This function depends on STATISTICS_WRAPPER
%
% TODO change cfg.frequency in all functions to cfg.foi or cfg.foilim

% Copyright (C) 2005-2006, Robert Oostenveld
%
% $Log: freqstatistics.m,v $
% Revision 1.23  2009/10/07 10:03:30  jansch
% temporary workaround to work with private copy statistics_wrapperJM, in order
% to develop some code. users whose (part of their) username contains 'jan'
% will run into problems, and have to uncomment lines 140, 143-145
%
% Revision 1.22  2009/04/08 15:57:08  roboos
% moved the handling of the output cfg (with all history details) from wrapper to main function
%
% Revision 1.21  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.20  2007/07/16 16:02:13  roboos
% fixed small bug for cfg.parameter
%
% Revision 1.19  2007/07/04 08:21:02  roboos
% remove cross/coherence spectrum in case cfg.parameter=powspctrm
%
% Revision 1.18  2007/06/13 06:19:40  roboos
% fixed bug in strcmp(cfg.parameter, ...), thanks to Vladimir
%
% Revision 1.17  2007/05/14 08:28:05  roboos
% changed handling of non-powspctrm input, also support cohspctrm
%
% Revision 1.16  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.15  2007/04/02 14:33:44  roboos
% disabled checkinput for the time being, since freq data can contain stat/zvalue/tvalue instead of powsptrm
%
% Revision 1.14  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.13  2007/03/27 15:20:50  erimar
% Updated help (replaced "p-value" by "significance probability").
%
% Revision 1.12  2007/01/17 13:23:35  roboos
% added bug report, no functional change
%
% Revision 1.11  2006/11/27 15:38:20  roboos
% implemented support for cfg.parameter, by locally renaming the field in the data structure
%
% Revision 1.10  2006/10/19 15:05:51  roboos
% updated documentation
%
% Revision 1.9  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.8  2006/07/12 09:18:17  roboos
% improved documentation
%
% Revision 1.7  2006/06/20 16:25:58  ingnie
% updated documentation
%
% Revision 1.6  2006/06/20 12:57:26  roboos
% updated documentation
%
% Revision 1.5  2006/06/13 14:48:09  ingnie
% updated documentation

fieldtripdefs

% check if the input data is valid for this function
for i=1:length(varargin)
  % FIXME at this moment (=2 April) this does not work, because the input might not always have a powspctrm o.i.d.
  % See email from Juriaan
  % varargin{i} = checkdata(varargin{i}, 'datatype', 'freq', 'feedback', 'no');
end

% the low-level data selection function does not know how to deal with other parameters, so work around it
if isfield(cfg, 'parameter') && strcmp(cfg.parameter, 'powspctrm')
  % use the power spectrum, this is the default 
  for i=1:length(varargin)
    if isfield(varargin{i}, 'crsspctrm'), varargin{i} = rmfield(varargin{i}, 'crsspctrm'); end % remove to avoid confusion
    if isfield(varargin{i}, 'cohspctrm'), varargin{i} = rmfield(varargin{i}, 'cohspctrm'); end % remove to avoid confusion
    if isfield(varargin{i}, 'labelcmb'),  varargin{i} = rmfield(varargin{i}, 'labelcmb');  end % remove to avoid confusion
  end
elseif isfield(cfg, 'parameter') && strcmp(cfg.parameter, 'crsspctrm')
  % use the cross spectrum, this might work as well (but has not been tested)
elseif isfield(cfg, 'parameter') && strcmp(cfg.parameter, 'cohspctrm')
  % for testing coherence on group level: 
  % rename cohspctrm->powspctrm and labelcmb->label
  for i=1:length(varargin)
    dat         = varargin{i}.cohspctrm;
    labcmb      = varargin{i}.labelcmb;
    for j=1:size(labcmb)
      lab{j,1} = sprintf('%s - %s', labcmb{j,1}, labcmb{j,2});
    end
    varargin{i} = rmsubfield(varargin{i}, 'cohspctrm');
    varargin{i} = rmsubfield(varargin{i}, 'labelcmb');
    varargin{i} = setsubfield(varargin{i}, 'powspctrm', dat);
    varargin{i} = setsubfield(varargin{i}, 'label',     lab);
  end
elseif isfield(cfg, 'parameter')
  % rename the desired parameter to powspctrm
  fprintf('renaming parameter ''%s'' into ''powspctrm''\n', cfg.parameter);
  for i=1:length(varargin)
    dat         = getsubfield(varargin{i}, cfg.parameter);
    varargin{i} = rmsubfield (varargin{i}, cfg.parameter);
    varargin{i} = setsubfield(varargin{i}, 'powspctrm', dat);
  end
end

[status,output] = system('whoami');
if isempty(strfind(output,'jan')),
  % call the general function
  [stat, cfg] = statistics_wrapper(cfg, varargin{:});
else
  % call the general function
  [stat, cfg] = statistics_wrapperJM(cfg, varargin{:});
end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: freqstatistics.m,v 1.23 2009/10/07 10:03:30 jansch Exp $';

% remember the configuration of the input data
cfg.previous = [];
for i=1:length(varargin)
  if isfield(varargin{i}, 'cfg')
    cfg.previous{i} = varargin{i}.cfg;
  else
    cfg.previous{i} = [];
  end
end

% remember the exact configuration details
stat.cfg = cfg;
