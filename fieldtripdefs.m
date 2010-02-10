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

% this is for Matlab version specific backward compatibility support
% the version specific path should only be added once in every session
persistent versionpath
persistent signalpath

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

if isempty(versionpath)
  % fieldtrip/compat contains version specific subdirectories to facilitate in backward compatibility support
  switch version('-release')
    case '13'
      % Version 6.5.1.199709 Release 13 (Service Pack 1)
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', 'R13');
    case '14'
      % Version 7.0.4.352 (R14) Service Pack 2
      % Version 7.1.0.183 (R14) Service Pack 3
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', 'R14');
    case '2006a'
      % Version 7.2.0.283 (R2006a)
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', '2006a');
    case '2006b'
      % Version 7.3
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', '2006b');
    case '2007a'
      % Version 7.4
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', '2007a');
    case '2007b'
      % Version 7.5
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', '2007b');
    case '2008a'
      % Version 7.6
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', '2008a');
    case '2008b'
      % Version 7.7
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', '2008b');
    case '2009a'
      % Version 7.8
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', '2009a');
    case '2009b'
      % Version 7.9
      versionpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', '2009b');
    otherwise
      % unknown release number, do nothing
      versionpath = 'unknown';
  end % switch
  if ~strcmp(versionpath, 'unknown') && isdir(versionpath) && ~isempty(strfind(path, versionpath))
    % only add the directory if it exists and was not added to the path before
    addpath(versionpath);
  end
end % if isempty(versionpath)


if isempty(signalpath)
  % test whether the signal processing toolbox is available
  if ~hastoolbox('signal')
    % add the fieldtrip/compat/signal directory to the path, which contains
    % some drop-in replacement code from the Octave project
    signalpath = fullfile(fileparts(which('fieldtripdefs')), 'compat', 'signal');
    addpath(signalpath);
  else
    % remember the location of the Mathworks signal processing toolbox
    signalpath = fileparts(which('butter'));
  end
end

