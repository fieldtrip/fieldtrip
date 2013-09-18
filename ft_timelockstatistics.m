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
%                    'montecarlo'    get Monte-Carlo estimates of the significance probabilities and/or critical values from the permutation distribution,
%                    'analytic'      get significance probabilities and/or critical values from the analytic reference distribution (typically, the sampling distribution under the null hypothesis),
%                    'stats'         use a parametric test from the Matlab statistics toolbox,
%                    'crossvalidate' use crossvalidation to compute predictive performance
%
% The other cfg options depend on the method that you select. You
% should read the help of the respective subfunction FT_STATISTICS_XXX
% for the corresponding configuration options and for a detailed
% explanation of each method.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_TIMELOCKANALYSIS, FT_TIMELOCKGRANDAVERAGE

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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar varargin

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'no');
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

% call the general function which does the actual work
[stat, cfg] = statistics_wrapper(cfg, varargin{:});

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous varargin
ft_postamble history stat
ft_postamble savevar stat
