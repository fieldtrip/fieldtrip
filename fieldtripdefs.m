function fieldtripdefs

% FIELDTRIPDEFS is called at the begin of all FieldTrip functions and
% contains some defaults and path settings
%
% Note that this should be a function and not a script, otherwise the
% hastoolbox function appears not be found in fieldtrip/private.

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the global defaults, the checkconfig function will copy these into the local configurations
global ft_default
if ~isfield(ft_default, 'trackconfig'), ft_default.trackconfig = 'off';    end % cleanup, report, off
if ~isfield(ft_default, 'checkconfig'), ft_default.checkconfig = 'loose';  end % pedantic, loose, silent
if ~isfield(ft_default, 'checksize'),   ft_default.checksize   = 1e5;      end % number in bytes, can be inf

if isempty(which('hastoolbox'))
  % the fieldtrip/public directory contains the hastoolbox function
  % which is required for the remainder of this script
  addpath(fullfile(fileparts(which('fieldtripdefs')), 'public'));
end

try
  % this directory contains the backward compatibility wrappers for the ft_xxx function name change
  hastoolbox('compat', 1, 1);
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
  % numerous functions depend on this module
  hastoolbox('fileio', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('forwinv', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('preproc', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('plotting', 1, 1);
end

try
  % this contains some examples for realtime processing
  hastoolbox('realtime', 1, 1);
end

