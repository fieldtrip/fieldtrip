function [stat] = ft_timelockstatistics(cfg, varargin)

% FT_TIMELOCKSTATISTICS  computes significance probabilities and/or critical values of a parametric statistical test
% or a non-parametric permutation test.
%
% Use as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
% where the input data is the result from either FT_TIMELOCKANALYSIS or
% FT_TIMELOCKGRANDAVERAGE.
%
% The configuration can contain the following options for data selection
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
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
% See also FT_TIMELOCKANALYSIS, FT_TIMELOCKGRANDAVERAGE

% Undocumented local options:
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%                     The data should be provided in a cell array
%   cfg.outputfile = one can specify output as file to save to disk


% This function depends on STATISTICS_WRAPPER

% Copyright (C) 2005-2006, Robert Oostenveld
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

% set the defaults
if ~isfield(cfg, 'inputfile'),    cfg.inputfile = [];          end
if ~isfield(cfg, 'outputfile'),   cfg.outputfile = [];         end

hasdata = nargin>1;
if ~isempty(cfg.inputfile) % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    for i=1:numel(cfg.inputfile)
      varargin{i} = loadvar(cfg.inputfile{i}, 'data'); % read datasets from array inputfile
    end
  end
end

% % check if the input data is valid for this function
% for i=1:length(varargin)
%   % FIXME at this moment (=10 May) this does not work, because the input might not always have an avg
%   % See freqstatistics
%   %varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'no');
% end

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
cfg.version.id = '$Id$';

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

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', stat); % use the variable name "data" in the output file
end

