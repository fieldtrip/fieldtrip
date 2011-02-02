function ft_defaults

% FT_DEFAULTS is called at the begin of all FieldTrip functions and
% contains some defaults and path settings
% (formerly known as fieldtripdefs.m)
%
% Note that this should be a function and not a script, otherwise the
% ft_hastoolbox function appears not be found in fieldtrip/private.

% Copyright (C) 2009, Robert Oostenveld
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
% $Id: ft_defaults.m 2017 2010-11-02 08:19:11Z jansch $

% set the global defaults, the ft_checkconfig function will copy these into the local configurations
global ft_default
if ~isfield(ft_default, 'trackconfig'), ft_default.trackconfig = 'off';    end % cleanup, report, off
if ~isfield(ft_default, 'checkconfig'), ft_default.checkconfig = 'loose';  end % pedantic, loose, silent
if ~isfield(ft_default, 'checksize'),   ft_default.checksize   = 1e5;      end % number in bytes, can be inf

% this is for Matlab version specific backward compatibility support
% the version specific path should only be added once in every session
persistent versionpath
persistent signalpath

% some people mess up their path settings with addpath(genpath(...))which
% results in different versions of SPM or other other toolboxes on the path
list = which('spm', '-all');
if length(list)>1
  [ws warned] = warning_once('multiple versions of SPM on your path will confuse FieldTrip');
  if warned % only throw the warning once
      for i=1:length(list)
          warning('one version of SPM is found here: %s', list{i});
      end
  end
end

if isempty(which('ft_hastoolbox'))
  % the fieldtrip/utilities directory contains the ft_hastoolbox function
  % which is required for the remainder of this script
  addpath(fullfile(fileparts(which('ft_defaults')), 'utilities'));
end

try
  % this directory contains the backward compatibility wrappers for the ft_xxx function name change
  ft_hastoolbox('compat', 3, 1); % not required
end

try 
  % this directory contains the backward compatibility wrappers for the fieldtrip/utilities functions
  ft_hastoolbox('utilities/compat', 3, 1);
end

try
  % this contains layouts and cortical meshes
  ft_hastoolbox('template', 1, 1);
end

try
  % this is used in statistics
  ft_hastoolbox('statfun', 1, 1);
end

try
  % this is used in definetrial
  ft_hastoolbox('trialfun', 1, 1);
end

try
  % this contains the low-level reading functions
  ft_hastoolbox('fileio', 1, 1);
  ft_hastoolbox('fileio/compat', 3, 1); % not required
end

try
  % this is for filtering time-series data
  ft_hastoolbox('preproc', 1, 1);
  ft_hastoolbox('preproc/compat', 3, 1); % not required
end

try
  % this contains forward models for the EEG and MEG volume conduction problem
  ft_hastoolbox('forward', 1, 1);
  ft_hastoolbox('forward/compat', 3, 1); % not required
end

try
  % numerous functions depend on this module
  ft_hastoolbox('inverse', 1, 1);
end

try
  % this contains intermediate-level plotting functions, e.g. multiplots and 3-d objects
  ft_hastoolbox('plotting', 1, 1);
  ft_hastoolbox('plotting/compat', 1, 1);
end

try 
  % this contains the functions to compute connecitivy metrics
  ft_hastoolbox('connectivity', 1,1);
end

try
  % this contains specific code and examples for realtime processing
  ft_hastoolbox('realtime', 3, 1);             % not required
  ft_hastoolbox('realtime/datasource', 3, 1);  % not required
end

try 
  % this contains intermediate-level functions for spectral analysis
  ft_hastoolbox('specest', 1, 1);
end
