function fieldtripdefs

% FIELDTRIPDEFS is called at the begin of all FieldTrip functions and
% contains some defaults and path settings
%
% Note that this should be a function and not a script, otherwise the
% hastoolbox function appears not be found in fieldtrip/private.

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
% $Id$

% set the global defaults, the checkconfig function will copy these into the local configurations
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
  warning('multiple versions of SPM on your path will confuse FieldTrip');
  for i=1:length(list)
    warning('spm is found here: %s', list{i});
  end
end

if isempty(which('hastoolbox'))
  % the fieldtrip/public directory contains the hastoolbox function
  % which is required for the remainder of this script
  addpath(fullfile(fileparts(which('fieldtripdefs')), 'public'));
end

try
  % this directory contains the backward compatibility wrappers for the ft_xxx function name change
  hastoolbox('compat', 3, 1); % not required
end

try
  % this contains layouts and cortical meshes
  hastoolbox('template', 1, 1);
end

try
  % this is used in statistics
  hastoolbox('statfun', 1, 1);
end

try
  % this is used in definetrial
  hastoolbox('trialfun', 1, 1);
end

try
  % this contains the low-level reading functions
  hastoolbox('fileio', 1, 1);
  hastoolbox('fileio/compat', 3, 1); % not required
end

try
  % this is for filtering time-series data
  hastoolbox('preproc', 1, 1);
  hastoolbox('preproc/compat', 3, 1); % not required
end

try
  % this contains forward models for the EEG and MEG volume conduction problem
  hastoolbox('forward', 1, 1);
  hastoolbox('forward/compat', 3, 1); % not required
end

try
  % numerous functions depend on this module
  hastoolbox('forwinv', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('inverse', 1, 1);
end

try
  % this contains intermediate-level plotting functions, e.g. multiplots and 3-d objects
  hastoolbox('plotting', 1, 1);
  hastoolbox('plotting/compat', 1, 1);
end

try 
  % this contains the functions to compute connecitivy metrics
  hastoolbox('connectivity', 1,1);
end

try
  % this contains specific code and examples for realtime processing
  hastoolbox('realtime', 3, 1);             % not required
  hastoolbox('realtime/datasource', 3, 1);  % not required
end

try 
  % this contains intermediate-level functions for spectral analysis
  hastoolbox('specest', 1, 1);
end
