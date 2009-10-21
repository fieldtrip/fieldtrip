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
% $Log: sourcestatistics.m,v $
% Revision 1.44  2009/10/07 10:06:21  jansch
% temporary workaround to work with private copy statistics_wrapperJM, in order to develop some code (of coures with the eventual benefit for all). users whose (part of their) username contains 'jan' may run into problems, and have to uncomment lines 157, 161-163
%
% Revision 1.43  2009/04/08 15:57:08  roboos
% moved the handling of the output cfg (with all history details) from wrapper to main function
%
% Revision 1.42  2008/12/05 14:47:05  ingnie
% updated help on cfg.roi
%
% Revision 1.41  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.40  2008/07/31 16:22:52  roboos
% added documentation pertaining to atlas ROIs and added check on input in case of ROI
%
% Revision 1.39  2007/05/30 07:08:08  roboos
% use checkdata instead of fixinside
%
% Revision 1.38  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.37  2007/04/03 10:04:20  roboos
% allow both source and volume data as input
%
% Revision 1.36  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.35  2006/10/19 15:05:51  roboos
% updated documentation
%
% Revision 1.34  2006/07/12 09:17:58  roboos
% improved documentation, fixed small bug in converting cfg option
%
% Revision 1.33  2006/07/05 10:21:56  roboos
% updaed documentation for consistency
%
% Revision 1.32  2006/06/13 14:48:10  ingnie
% updated documentation
%
% Revision 1.31  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.30  2006/03/30 12:24:33  roboos
% Implemented private/fixinside, which facilitates consistent
% handling of source/volume data. Improved documentation. Fixed some
% bugs related to inconsistent handling of ROIs (i.e. inside/outside)
%
% Revision 1.29  2005/10/05 11:06:30  roboos
% documentation change
%
% Revision 1.28  2005/05/17 17:50:39  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.27  2005/05/12 07:20:06  roboos
% fixed backward compatibility approach/method translation for randcluster
%
% Revision 1.26  2005/05/09 14:20:02  roboos
% made the check on cfg.method optional, so that it also works when only the approach is specified
%
% Revision 1.25  2005/04/07 17:12:21  roboos
% improved the general help documentation, no code changes
%
% Revision 1.24  2005/04/06 14:29:41  roboos
% added copyrights and a log message placeholder
%

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
