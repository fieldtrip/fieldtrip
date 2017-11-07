function [argout, optout] = fexec(argin, optin)

% FEXEC is the low-level function that executes the job on the engine or
% slave. It also tries to change the path and pwd to those on the master
% and it catches and deals with any errors in the code that is executed.
%
% This function should not be called directly.
%
% See also PEEREXEC, QSUBEXEC, ENGEXEC

% -----------------------------------------------------------------------
% Copyright (C) 2011-2012, Robert Oostenveld
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/
% -----------------------------------------------------------------------

% keep track of the time
stopwatch = tic;

% these variables will be used when an error is caught
usediary  = false;
diaryfile = '';

argout = {};
optout = {};

% clear the previous warning and error messages
lastwarn('');
lasterr('');

% these will be determined later on, but are set here to empty for better error handling
masterid = [];
timallow = [];
memallow = [];

% there are many reasons why the execution may fail, hence the elaborate try-catch
try
  
  if ~iscell(argin)
    error('input argument should be a cell-array');
  end
  
  if ~ischar(argin{1}) && ~isa(argin{1}, 'function_handle')
    error('input argument #1 should be a string or a function handle');
  end
  
  fname = argin{1};
  argin = argin(2:end);
  
  if ~iscell(optin)
    error('input options should be a cell-array');
  end
  
  % check whether a diary file should be created
  usediary = ft_getopt(optin, 'diary');
  usediary = any(strcmp(usediary, {'always', 'warning', 'error'}));
  
  % check whether a watchdog should be set
  % this only applies to the peer distributed computing system
  masterid = ft_getopt(optin, 'masterid');
  timallow = ft_getopt(optin, 'timallow');
  memallow = []; % ft_getopt(optin, 'memallow');
  if ~isempty(masterid) || ~isempty(timallow) || ~isempty(memallow)
    watchdog(masterid, timallow, memallow);
  end
  
  % try setting the same path directory
  option_path = ft_getopt(optin, 'path');
  setcustompath(option_path);
  
  % try changing to the same working directory
  option_pwd = ft_getopt(optin, 'pwd');
  setcustompwd(option_pwd);
  
  % try assigning the same global variables
  option_global = ft_getopt(optin, 'global');
  setglobal(option_global);
  
  % seed the random number generator
  option_randomseed = ft_getopt(optin, 'randomseed');
  if ~isempty(option_randomseed)
    if ft_platform_supports('RandStream.setGlobalStream')
      % version 2012a gives a warning that RandStream.setDefaultStream will be removed in the future
      % and that RandStream.setGlobalStream should be used instead
      s = RandStream('mcg16807', 'Seed', option_randomseed);
      RandStream.setGlobalStream(s);
    elseif ft_platform_supports('RandStream.setDefaultStream')
      % this is according to http://www.mathworks.com/help/techdoc/math/bsn94u0-1.html
      % and is needed to avoid a warning about Using 'seed' to set RAND's internal state causes RAND, RANDI, and RANDN to use legacy random number generators.
      s = RandStream('mcg16807', 'Seed', option_randomseed);
      RandStream.setDefaultStream(s);
    else
      % in older Matlab versions, and in GNU Octave, it works like this
      rand ('seed', option_randomseed);
      randn('seed', option_randomseed);
    end
  end
  
  % there are potentially errors to catch from the which() function
  if ischar(fname) && isempty(which(fname))
    error('Not a valid M-file (%s).', fname);
  end
  
  % this controls the amount of data that has to be sent back, furthermore the
  % function internal operations and its output can depend on the number of
  % output arguments
  numargout = ft_getopt(optin, 'nargout');
  
  if isempty(numargout)
    % it can be difficult to determine the number of output arguments
    try
      if (isequal(fname, 'cellfun') || isequal(fname, @cellfun))
        numargout = nargout(argin{1});
      else
        numargout = nargout(fname);
      end
    catch
      % the "catch me" syntax is broken on MATLAB74, this fixes it
      nargout_error = lasterror;
      if strcmp(nargout_error.identifier, 'MATLAB:narginout:doesNotApply')
        % e.g. in case of nargin('plus')
        numargout = 1;
      else
        rethrow(nargout_error);
      end
    end
  end % determine number of output arguments
  
  if numargout<0
    % the nargout function returns -1 in case of a variable number of output arguments
    numargout = 1;
  end
  
  
  % start measuring the time and memory requirements
  memprofile on
  timused = toc(stopwatch);
  
  if usediary
    % switch on the diary
    diaryfile = tempname;
    diary(diaryfile);
  end
  
  % evaluate the function and get the output arguments
  argout = cell(1, numargout);
  if numargout>0
    [argout{:}] = feval(fname, argin{:});
  else
    feval(fname, argin{:});
  end
  
  if usediary && exist(diaryfile, 'file')
    % close the diary and read the contents
    diary off
    fid = fopen(diaryfile, 'r');
    diarystring = fread(fid, [1, inf], 'char=>char');
    fclose(fid);
  else
    % return an empty diary
    diarystring = [];
  end
  
  % determine the time and memory requirements
  timused = toc(stopwatch) - timused;
  memstat = memprofile('info');
  memprofile off
  memprofile clear
  
  % determine the maximum amount of memory that was used during the function evaluation
  memused = max([memstat.mem]) - min([memstat.mem]);
  
  % Note that the estimated memory is inaccurate, because of
  % the dynamic memory management of Matlab and the garbage
  % collector. Especially on small jobs, the reported memory
  % use does not replect the size of the variables involved in
  % the computation. Matlab is able to squeeze these small jobs
  % in some left-over memory fragment that was not yet deallocated.
  % Larger memory jobs return more reliable measurements.
  
  fprintf('executing job took %f seconds and %d bytes\n', timused, memused);
  
  % collect the output options
  optout = {'timused', timused, 'memused', memused, 'lastwarn', lastwarn, 'lasterr', '', 'diary', diarystring, 'release', version('-release'), 'pwd', pwd, 'path', path, 'hostname', getenv('HOSTNAME')};
  
catch
  % the "catch me" syntax is broken on MATLAB74, this fixes it
  feval_error = lasterror;
  
  if usediary && exist(diaryfile, 'file')
    % close the diary and read the contents
    diary off
    fid = fopen(diaryfile, 'r');
    diarystring = fread(fid, [1, inf], 'char=>char');
    fclose(fid);
  else
    % return an empty diary
    diarystring = [];
  end
  
  % the output options will include the error
  % note that the error cannot be sent as object, but has to be sent as struct
  optout = {'lastwarn', lastwarn, 'lasterr', struct(feval_error), 'diary', diarystring, 'release', version('-release'), 'pwd', pwd, 'path', path};
  
  % an error was detected while executing the job
  warning('an error was detected during job execution');
  
  % ensure that the memory profiler is switched off
  memprofile off
  memprofile clear
end % try-catch

% revert to the original path
% if ~isempty(option_path)
%   path(orig_path);
% end

% revert to the original working directory
% if ~isempty(option_pwd)
%   cd(orig_pwd);
% end

% clear the function and any persistent variables in it
if ischar(fname)
  clear(fname);
end

% close all files and figures
fclose all;
close all force;
close all hidden;

% clear any global variables
clear global

% clear the optional watchdog, which is loaded into memory as a mex file
if ~isempty(masterid) || ~isempty(timallow) || ~isempty(memallow)
  watchdog(0,0,0); % this is required to unlock it from memory
end

