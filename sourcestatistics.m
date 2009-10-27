function [stat] = sourcestatistics(cfg, varargin)

% SOURCESTATISTICS computes the probability for a given null-hypothesis using
% a parametric statistical test or using a non-parametric randomization test.
% 
% Use as
%   [stat] = sourcestatistics(cfg, source1, source2, ...)
% where the input data is the result from SOURCEANALYSIS, SOURCEDESCRIPTIVES
% or SOURCEGRANDAVERAGE.  The source structures should be spatially alligned
% to each other and should have the same positions for the source grid.
%
% The configuration should contain the following option for data selection
%   cfg.parameter  = string, describing the functional data to be processed, e.g. 'pow', 'nai' or 'coh'
%
% Furthermore, the configuration should contain:
%   cfg.method       = different methods for calculating the probability of the null-hypothesis,
%                    'montecarlo'    uses a non-parametric randomization test to get a Monte-Carlo estimate of the probability,
%                    'analytic'      uses a parametric test that results in analytic probability,
%                    'glm'           uses a general linear model approach,
%                    'stats'         uses a parametric test from the Matlab statistics toolbox,
%                    'parametric'    uses the Matlab statistics toolbox (very similar to 'stats'),
%                    'randomization' uses randomization of the data prior to source reconstruction,
%                    'randcluster'   uses randomization of the data prior to source reconstruction 
%                                    in combination with spatial clusters.
%
% You can restrict the statistical analysis to regions of interest (ROIs)
% or to the average value inside ROIs using the following options:
%   cfg.atlas        = filename of the atlas
%   cfg.roi          = string or cell of strings, region(s) of interest from anatomical atlas
%   cfg.avgoverroi   = 'yes' or 'no' (default = 'no')
%   cfg.hemisphere   = 'left', 'right', 'both', 'combined', specifying this is
%                      required when averaging over regions
%   cfg.inputcoord   = 'mni' or 'tal', the coordinate system in which your source 
%                      reconstruction is expressed
%
% The other cfg options depend on the method that you select. You
% should read the help of the respective subfunction STATISTICS_XXX
% for the corresponding configuration options and for a detailed
% explanation of each method.
%
%
% See also SOURCEANALYSIS, SOURCEDESCRIPTIVES, SOURCEGRANDAVERAGE

% Undocumented local options:
%   cfg.statistic

% Copyright (C) 2005-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

% this wrapper should be compatible with the already existing statistical
% functions that only work for source input data

% check if the input data is valid for this function
for i=1:length(varargin)
  if isfield(cfg, 'roi') && ~isempty(cfg.roi)
    varargin{i} = checkdata(varargin{i}, 'datatype', 'source', 'feedback', 'no', 'inside', 'index');
  else
    varargin{i} = checkdata(varargin{i}, 'datatype', {'source', 'volume'}, 'feedback', 'no', 'inside', 'index');
  end
end

if isfield(cfg, 'method')
  % call the appropriate subfunction
  if (strcmp(cfg.method, 'zero-baseline') || ...
      strcmp(cfg.method, 'nai')           || ...
      strcmp(cfg.method, 'pseudo-t')      || ...
      strcmp(cfg.method, 'difference')    || ...
      strcmp(cfg.method, 'anova1')        || ...
      strcmp(cfg.method, 'kruskalwallis'))
    % these are all statistical methods that are implemented in the old SOURCESTATISTICS_PARAMETRIC subfunction
    cfg.statistic = cfg.method;
    cfg.method    = 'parametric';
  elseif strcmp(cfg.method, 'randomization')
    cfg.method = 'randomization';
  elseif strcmp(cfg.method, 'randcluster')
    cfg.method = 'randcluster';
  end
end

if strcmp(cfg.method, 'parametric')
  % use the source-specific statistical subfunction
  stat = sourcestatistics_parametric(cfg, varargin{:});
elseif strcmp(cfg.method, 'randomization')
  % use the source-specific statistical subfunction
  stat = sourcestatistics_randomization(cfg, varargin{:});
elseif strcmp(cfg.method, 'randcluster')
  % use the source-specific statistical subfunction
  stat = sourcestatistics_randcluster(cfg, varargin{:});
else
  [status,output] = system('whoami');
  if isempty(strfind(output,'jan')),
    % use the data-indepentend statistical wrapper function
    % this will collect the data and subsequently call STATISTICS_XXX
    [stat, cfg] = statistics_wrapper(cfg, varargin{:});
  else
    [stat, cfg] = statistics_wrapperJM(cfg, varargin{:});
  end
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
cfg.version.id = '$Id: sourcestatistics.m,v 1.44 2009/10/07 10:06:21 jansch Exp $';

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
