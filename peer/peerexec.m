function [argout, optout] = peerexec(argin, optin)

% PEEREXEC is the low-level function that executes the job on the
% slave. It also tries to change the path and pwd to those on the
% master and it catches and deals with any errors in the code that
% is executed.
%
% This function should not be called directly.
%
% See also PEERSLAVE, PEERMASTER, PEERCELLFUN

% -----------------------------------------------------------------------
% Copyright (C) 2010, Robert Oostenveld
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

% there are many reasons why the execution may fail, hence the elaborate try-catch
try

  if ~iscell(argin)
    error('input argument should be a cell-array');
  end

  if ~ischar(argin{1})
    error('input argument #1 should be a string');
  end

  fname = argin{1};
  argin = argin(2:end);

  if ~iscell(optin)
    error('input options should be a cell-array');
  end

  % check whether a kill switch should be set
  masterid = keyval('masterid', optin);
  timallow = keyval('timallow', optin);
  if ~isempty(masterid) || ~isempty(timallow)
    killswitch(masterid, time+timallow);
  end

  % check whether a diary file should be created
  usediary = keyval('diary', optin);
  usediary = any(strcmp(usediary, {'always', 'warning', 'error'}));

  % try setting the same path directory
  option_path = keyval('path', optin);
  setcustompath(option_path);

  % try changing to the same working directory
  option_pwd = keyval('pwd', optin);
  setcustompwd(option_pwd);

  % there are potentially errors to catch from the which() function
  if isempty(which(fname))
    error('Not a valid M-file (%s).', fname);
  end

  % it can be difficult to determine the number of output arguments
  try
    numargout = nargout(fname);
  catch nargout_error
    if strcmp(nargout_error.identifier, 'MATLAB:narginout:doesNotApply')
      % e.g. in case of nargin('plus')
      numargout = 1;
    else
      rethrow(nargout_error);
    end
  end

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
  argout  = cell(1, numargout);
  [argout{:}] = feval(fname, argin{:});

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
  optout = {'timused', timused, 'memused', memused, 'lastwarn', lastwarn, 'lasterr', '', 'diary', diarystring, 'release', version('-release')};

catch feval_error

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
  optout = {'lastwarn', lastwarn, 'lasterr', struct(feval_error), 'diary', diarystring};

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
clear(fname);

% close all files and figures
fclose all;
close all force;
close all hidden;

% clear any global variables
clear global

% clear the optional kill switch, which is loaded into memory as a mex file
clear killswitch

% clear the previous warning and error messages
lastwarn('');
lasterr('');

