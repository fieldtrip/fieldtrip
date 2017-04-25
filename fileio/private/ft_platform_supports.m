function tf = ft_platform_supports(what,varargin)

% FT_PLATFORM_SUPPORTS returns a boolean indicating whether the current platform
% supports a specific capability
%
% Use as
%   status = ft_platform_supports(what)
% or
%   status = ft_platform_supports('matlabversion', min_version, max_version)
%
% The following values are allowed for the 'what' parameter, which means means that
% the specific feature explained on the right is supported:
%
%   'which-all'                     which(...,'all')
%   'exists-in-private-directory'   exists(...) will look in the /private subdirectory to see if a file exists
%   'onCleanup'                     onCleanup(...)
%   'alim'                          alim(...)
%   'int32_logical_operations'      bitand(a,b) with a, b of type int32
%   'graphics_objects'              graphics sysem is object-oriented
%   'libmx_c_interface'             libmx is supported through mex in the C-language (recent MATLAB versions only support C++)
%   'stats'                         all statistical functions in FieldTrip's external/stats directory
%   'program_invocation_name'       program_invocation_name() (GNU Octave)
%   'singleCompThread'              start MATLAB with -singleCompThread
%   'nosplash'                      start MATLAB with -nosplash
%   'nodisplay'                     start MATLAB with -nodisplay
%   'nojvm'                         start MATLAB with -nojvm
%   'no-gui'                        start GNU Octave with --no-gui
%   'RandStream.setGlobalStream'    RandStream.setGlobalStream(...)
%   'RandStream.setDefaultStream'   RandStream.setDefaultStream(...)
%   'rng'                           rng(...)
%   'rand-state'                    rand('state')
%   'urlread-timeout'               urlread(..., 'Timeout', t)
%   'griddata-vector-input'         griddata(...,...,...,a,b) with a and b vectors
%   'griddata-v4'                   griddata(...,...,...,...,...,'v4') with v4 interpolation support
%   'uimenu'                        uimenu(...)
%   'weboptions'                    weboptions(...)
%   'parula'                        parula(...)
%
% See also FT_VERSION, VERSION, VER, VERLESSTHAN

if ~ischar(what)
  error('first argument must be a string');
end

switch what
  case 'matlabversion'
    tf = is_matlab() && matlabversion(varargin{:});
    
  case 'exists-in-private-directory'
    tf = is_matlab();
    
  case 'which-all'
    tf = is_matlab();
    
  case 'onCleanup'
    tf = is_octave() || matlabversion(7.8, Inf);
    
  case 'alim'
    tf = is_matlab();
    
  case 'int32_logical_operations'
    % earlier version of MATLAB don't support bitand (and similar)
    % operations on int32
    tf = is_octave() || ~matlabversion(-Inf, '2012a');
    
  case 'graphics_objects'
    % introduced in MATLAB 2014b, graphics is handled through objects;
    % previous versions use numeric handles
    tf = is_matlab() && matlabversion('2014b', Inf);
    
  case 'libmx_c_interface'
    % removed after 2013b
    tf = is_matlab() && matlabversion(-Inf, '2013b');
    
  case 'stats'
    root_dir = fileparts(which('ft_defaults'));
    external_stats_dir = fullfile(root_dir, 'external', 'stats');
    
    % these files are only used by other functions in the external/stats directory
    exclude_mfiles = {
      'common_size.m'
      'iscomplex.m'
      'lgamma.m'
      };
    
    tf = has_all_functions_in_dir(external_stats_dir, exclude_mfiles);
    
  case 'program_invocation_name'
    % Octave supports program_invocation_name, which returns the path
    % of the binary that was run to start Octave
    tf = is_octave();
    
  case 'singleCompThread'
    tf = is_matlab() && matlabversion(7.8, Inf);
    
  case {'nosplash', 'nodisplay', 'nojvm'}
    % Only on MATLAB
    tf = is_matlab();
    
  case 'no-gui'
    % Only on Octave
    tf = is_octave();
    
  case 'RandStream.setDefaultStream'
    tf = is_matlab() && matlabversion('2008b', '2011b');
    
  case 'RandStream.setGlobalStream'
    tf = is_matlab() && matlabversion('2012a', Inf);
    
  case 'randomized_PRNG_on_startup'
    tf = is_octave() || ~matlabversion(-Inf, '7.3');
    
  case 'rng'
    % recent MATLAB versions
    tf = is_matlab() && matlabversion('7.12', Inf);
    
  case 'rand-state'
    % GNU Octave
    tf = is_octave();
    
  case 'urlread-timeout'
    tf = is_matlab() && matlabversion('2012b', Inf);
    
  case 'griddata-vector-input'
    tf = is_matlab();
    
  case 'griddata-v4'
    tf = is_matlab() && matlabversion('2009a', Inf);
    
  case 'uimenu'
    tf = is_matlab();
    
  case {'weboptions', 'webread', 'websave'}
    tf = is_matlab() && matlabversion('2014b', Inf);
    
  case 'webwrite'
    tf = is_matlab() && matlabversion('2015a', Inf);
    
  case 'boundary'
    tf = is_matlab() && matlabversion('2014b', Inf);
    
  case 'parula'
    tf = is_matlab() && matlabversion('2014b', Inf);
    
  otherwise
    error('unsupported value for first argument: %s', what);
    
end % switch

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = is_matlab()
tf = ~is_octave();
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = is_octave()
persistent cached_tf

if isempty(cached_tf)
  cached_tf = logical(exist('OCTAVE_VERSION', 'builtin'));
end

tf = cached_tf;
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = has_all_functions_in_dir(in_dir, exclude_mfiles)
% returns true if all functions in in_dir are already provided by the platform
m_files = dir(fullfile(in_dir,'*.m'));
n = numel(m_files);

for k = 1:n
  m_filename = m_files(k).name;
  if isempty(which(m_filename)) && isempty(strmatch(m_filename,exclude_mfiles))
    tf = false;
    return;
  end
end

tf = true;

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inInterval] = matlabversion(min, max)

% MATLABVERSION checks if the current MATLAB version is within the interval
% specified by min and max.
%
% Use as
%   if matlabversion(7.0, 7.9)
%      % do something
%   end
%
% Both strings and numbers, as well as infinities, are supported, eg.:
%   matlabversion(7.1, 7.9)         % is version between 7.1 and 7.9?
%   matlabversion(6, '7.10')        % is version between 6 and 7.10? (note: '7.10', not 7.10)
%   matlabversion(-Inf, 7.6)        % is version <= 7.6?
%   matlabversion('2009b')          % exactly 2009b
%   matlabversion('2008b', '2010a') % between two versions
%   matlabversion('2008b', Inf)     % from a version onwards
% etc.
%
% See also VERSION, VER, VERLESSTHAN

% Copyright (C) 2006, Robert Oostenveld
% Copyright (C) 2010, Eelke Spaak
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
% FieldTrip is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% FieldTrip is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% the version does not change, making it persistent speeds up the subsequent calls
persistent curVer

if nargin<2
  max = min;
end

if isempty(curVer)
  curVer = version();
end

if ((ischar(min) && isempty(str2num(min))) || (ischar(max) && isempty(str2num(max))))
  % perform comparison with respect to release string
  
  ind = strfind(curVer, '(R');
  [year, ab] = parseMatlabRelease(curVer((ind + 2):(numel(curVer) - 1)));
  
  [minY, minAb] = parseMatlabRelease(min);
  [maxY, maxAb] = parseMatlabRelease(max);
  
  inInterval = orderedComparison(minY, minAb, maxY, maxAb, year, ab);
  
else
  % perform comparison with respect to version number
  [major, minor] = parseMatlabVersion(curVer);
  [minMajor, minMinor] = parseMatlabVersion(min);
  [maxMajor, maxMinor] = parseMatlabVersion(max);
  
  inInterval = orderedComparison(minMajor, minMinor, maxMajor, maxMinor, major, minor);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [year, ab] = parseMatlabRelease(str)
if (str == Inf)
  year = Inf; ab = Inf;
elseif (str == -Inf)
  year = -Inf; ab = -Inf;
else
  year = str2num(str(1:4));
  ab = str(5);
end
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [major, minor] = parseMatlabVersion(ver)
if (ver == Inf)
  major = Inf; minor = Inf;
elseif (ver == -Inf)
  major = -Inf; minor = -Inf;
elseif (isnumeric(ver))
  major = floor(ver);
  minor = int8((ver - floor(ver)) * 10);
else % ver is string (e.g. '7.10'), parse accordingly
  [major, rest] = strtok(ver, '.');
  major = str2num(major);
  minor = str2num(strtok(rest, '.'));
end
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inInterval = orderedComparison(lowerA, lowerB, upperA, upperB, testA, testB)
% checks if testA is in interval (lowerA,upperA); if at edges, checks if testB is in interval (lowerB,upperB).
if (testA < lowerA || testA > upperA)
  inInterval = false;
else
  inInterval = true;
  if (testA == lowerA)
    inInterval = inInterval && (testB >= lowerB);
  end
  
  if (testA == upperA)
    inInterval = inInterval && (testB <= upperB);
  end
end
end % function
