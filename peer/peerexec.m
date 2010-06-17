function [argout, options] = peerexec(argin, options)

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

usediary = false;
% keep track of the time and number of jobs
stopwatch = tic;
prevtime  = toc(stopwatch);
idlestart = toc(stopwatch);
jobnum    = 0;

    
try
  % there are many reasons why the execution may fail, hence the elaborate try-catch
  
  if ~iscell(argin)
    error('input argument should be a cell-array');
  end
  
  if ~ischar(argin{1})
    error('input argument #1 should be a string');
  end
  
  fname = argin{1};
  argin = argin(2:end);
  
  if ~iscell(options)
    error('input options should be a cell-array');
  end
  
  % try setting the same path directory
  option_path = keyval('path', options);
  if ~isempty(option_path)
    path(option_path, path);
  end
  
  % try changing to the same working directory
  option_pwd = keyval('pwd', options);
  if ~isempty(option_pwd)
    try
      cd(option_pwd);
    catch cd_error
      % don't throw an error, just give a warning (and hope for the best...)
      warning(cd_error.message);
    end
  end
  
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
  
  if usediary
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
  
  fprintf('executing job %d took %f seconds and %d bytes\n', jobnum, timused, memused);
  
  % collect the output options
  options = {'timused', timused, 'memused', memused, 'lastwarn', lastwarn, 'lasterr', '', 'diary', diarystring, 'release', version('-release')};
  
catch feval_error
  argout  = {};
  
  if usediary
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
  options = {'lastwarn', lastwarn, 'lasterr', feval_error, 'diary', diarystring};
  % an error was detected while executing the job
  warning('an error was detected during job execution');
  % ensure that the memory profiler is switched off
  memprofile off
  memprofile clear
end % try-catch

