function [stat] = timelockstatistics(cfg, varargin)

% TIMELOCKSTATISTICS  computes significance probabilities and/or critical values of a parametric statistical test
% or a non-parametric permutation test.
%
% Use as
%   [stat] = timelockstatistics(cfg, timelock1, timelock2, ...)
% where the input data is the result from either TIMELOCKANALYSIS or
% TIMELOCKGRANDAVERAGE.
%
% The configuration can contain the following options for data selection
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see CHANNELSELECTION for details
%   cfg.latency     = [begin end] in seconds or 'all' (default = 'all')
%   cfg.avgoverchan = 'yes' or 'no'                   (default = 'no')
%   cfg.avgovertime = 'yes' or 'no'                   (default = 'no')
%   cfg.parameter   = string                          (default = 'trial' or 'avg')
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
% See also TIMELOCKANALYSIS, TIMELOCKGRANDAVERAGE

% This function depends on STATISTICS_WRAPPER

% Copyright (C) 2005-2006, Robert Oostenveld
%
% $Log: timelockstatistics.m,v $
% Revision 1.25  2009/04/08 15:57:08  roboos
% moved the handling of the output cfg (with all history details) from wrapper to main function
%
% Revision 1.24  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.23  2007/05/10 10:18:39  ingnie
% disabled checkinput for the time being, since timelock data can contain stat/zvalue/tvalue instead of an avg
%
% Revision 1.22  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.21  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.20  2007/03/27 15:19:14  erimar
% Updated help (replaced "p-value" by "significance probability".
%
% Revision 1.19  2006/11/27 15:38:20  roboos
% implemented support for cfg.parameter, by locally renaming the field in the data structure
%
% Revision 1.18  2006/10/19 15:05:51  roboos
% updated documentation
%
% Revision 1.17  2006/10/04 07:10:08  roboos
% updated documentation
%
% Revision 1.16  2006/07/12 09:18:17  roboos
% improved documentation
%
% Revision 1.15  2006/06/20 12:57:51  roboos
% updated documentation, removed support for Jens' old options
%
% Revision 1.14  2006/06/13 14:48:12  ingnie
% updated documentation

fieldtripdefs

% check if the input data is valid for this function
for i=1:length(varargin)
  % FIXME at this moment (=10 May) this does not work, because the input might not always have an avg
  % See freqstatistics
  %varargin{i} = checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'no');
end

% the low-level data selection function does not know how to deal with other parameters, so work around it
if isfield(cfg, 'parameter')
  if strcmp(cfg.parameter, 'trial') || strcmp(cfg.parameter, 'individual')
    % this is dealt with correctly in the low-level code, even if an average is present
  elseif strcmp(cfg.parameter, 'avg')
    % this is only dealt with in the low-level code if no single-trial/individual data is present
    for i=1:length(varargin)
      if isfield(varargin{i}, 'trial')
        varargin{i} = rmfield(varargin{i}, 'trial');
      end
      if isfield(varargin{i}, 'individual')
        varargin{i} = rmfield(varargin{i}, 'individual');
      end
    end
  else
    % rename the parameter of interest into 'avg'
    fprintf('renaming parameter ''%s'' into ''avg''\n', cfg.parameter);
    for i=1:length(varargin)
      dat         = getsubfield(varargin{i}, cfg.parameter);
      varargin{i} = rmsubfield (varargin{i}, cfg.parameter);
      varargin{i} = setsubfield(varargin{i}, 'avg', dat);
      if isfield(varargin{i}, 'trial')
        varargin{i} = rmfield(varargin{i}, 'trial');
      end
      if isfield(varargin{i}, 'individual')
        varargin{i} = rmfield(varargin{i}, 'individual');
      end
    end
  end
end

% call the general function
[stat, cfg] = statistics_wrapper(cfg, varargin{:});

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: timelockstatistics.m,v 1.25 2009/04/08 15:57:08 roboos Exp $';

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
